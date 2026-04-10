"""
HTTP endpoints for the test workflow.

The data flow is:

    GET  /api/tests                  → list all stored test records
    POST /api/tests                  → save a single record
    GET  /api/tests/load             → fetch one record (falls back to disk)
    GET  /api/tests/file-paths       → resolve where a test would be written
    POST /api/tests/write            → write a test .cc file to disk
    POST /api/tests/compile-run      → write + compile + run, stream all output
    POST /api/tests/git-commit       → git add + commit a test file
    POST /api/tests/validate-fix     → compile/run with auto-fix retries
    POST /api/tests/generate-all     → batch: oracle → generate → fix for every fn

    POST /api/generate-test          → single function: oracle → stream test code

The interesting endpoints (`/api/generate-test`, `/api/tests/validate-fix`,
`/api/tests/generate-all`) all delegate to the helpers in `pipeline.py`, so
the actual oracle/llm/compile/fix logic lives in one place.
"""
import asyncio
import json
import os
from pathlib import Path

from fastapi import APIRouter, Query
from fastapi.responses import StreamingResponse

import report
from compiler import make_compile_cmd
from config import MAX_TEST_RETRIES, PROBE_PDB_PATH, args
from models import (
    CompileRunRequest,
    GenerateAllRequest,
    GenerateTestRequest,
    GitCommitRequest,
    TestRecord,
    ValidateFixRequest,
    WriteTestRequest,
)
from ollama import call_ollama
from pipeline import (
    compile_run_fix_loop,
    oracle_output_text,
    oracle_probe_source,
    run_oracle_for_function,
    stream_test_generation,
)
from prompts import build_generate_test, parse_both_sections, strip_fences
from storage import load_tests, save_tests
from test_utils import get_test_file_path, wrap_test_content

router = APIRouter()


def _ndjson(event: dict) -> str:
    """Encode an event as one NDJSON line."""
    return json.dumps(event) + "\n"


# ── Basic CRUD ────────────────────────────────────────────────────────────────

@router.get("/api/tests")
def get_tests():
    return load_tests()


@router.post("/api/tests")
def save_test(record: TestRecord):
    data = load_tests()
    data[record.key] = {
        "mmdb_test":  record.mmdb_test,
        "gemmi_test": record.gemmi_test,
        "notes":      record.notes,
        "status":     record.status,
    }
    save_tests(data)
    return {"ok": True, "key": record.key}


@router.get("/api/tests/load")
def load_test_for_function(
    rel_source_path: str = Query(...),
    fn_name: str = Query(...),
    fn_line: int = Query(...),
):
    """Return the stored record for a function, falling back to .cc files on disk."""
    key = f"{rel_source_path}::{fn_name}:{fn_line}"
    record = load_tests().get(key, {})

    mmdb_content  = record.get("mmdb_test", "")
    gemmi_content = record.get("gemmi_test", "")

    if not mmdb_content:
        p = get_test_file_path(rel_source_path, fn_name, "mmdb")
        if p.exists():
            mmdb_content = p.read_text(errors="replace")

    if not gemmi_content:
        p = get_test_file_path(rel_source_path, fn_name, "gemmi")
        if p.exists():
            gemmi_content = p.read_text(errors="replace")

    return {
        "mmdb_test":  mmdb_content,
        "gemmi_test": gemmi_content,
        "notes":      record.get("notes", ""),
        "status":     record.get("status", "draft"),
    }


# ── Generate (single function, streaming) ────────────────────────────────────

@router.post("/api/generate-test")
async def generate_test(req: GenerateTestRequest):
    """Generate a test for one function.

    For MMDB-involving targets we first run an oracle probe to capture the
    function's real output on a real PDB file, then feed those measured
    values into the test-generation prompt as expected values. The endpoint
    streams NDJSON so the frontend can show probe progress AND the LLM's
    thinking + token-by-token test code.

    Event types:
        oracle_start | oracle_event | oracle_done | oracle_skip
        test_start   | test_thinking | test_chunk | test_done
        error
    """
    run_oracle = req.target in ("mmdb", "both")

    async def stream_events():
        oracle_result = None

        # ── Phase 1: oracle probe (optional) ─────────────────────────────
        if run_oracle:
            async for evt in run_oracle_for_function(
                model=req.model,
                function_name=req.function_name,
                source_code=req.source_code,
                mmdb_symbols=req.mmdb_symbols,
                additional_instructions=req.additional_instructions,
                rel_source_path=req.rel_source_path,
            ):
                if evt["type"] == "pipeline_result":
                    oracle_result = evt["result"]
                    continue
                yield _ndjson(evt)

        # ── Phase 2: stream the generated test code ──────────────────────
        async for evt in stream_test_generation(
            model=req.model,
            function_name=req.function_name,
            source_code=req.source_code,
            mmdb_symbols=req.mmdb_symbols,
            target=req.target,
            additional_instructions=req.additional_instructions,
            oracle_output=oracle_output_text(oracle_result),
            rel_source_path=req.rel_source_path,
            probe_source=oracle_probe_source(oracle_result),
            probe_pdb_path=PROBE_PDB_PATH,
        ):
            # `test_done` carries the full assembled text — drop it from the
            # wire to keep the existing client contract (it already accumulated
            # the chunks itself).
            if evt["type"] == "test_done":
                yield _ndjson({"type": "test_done"})
                continue
            yield _ndjson(evt)

    return StreamingResponse(stream_events(), media_type="application/x-ndjson")


# ── File paths / write / compile-run ─────────────────────────────────────────

@router.get("/api/tests/file-paths")
def get_test_file_paths(rel_source_path: str = Query(...), fn_name: str = Query(...)):
    """Return where MMDB and Gemmi test files would be written."""
    mmdb_path  = get_test_file_path(rel_source_path, fn_name, "mmdb")
    gemmi_path = get_test_file_path(rel_source_path, fn_name, "gemmi")
    return {
        "mmdb":         str(mmdb_path.relative_to(report._coot_root)),
        "gemmi":        str(gemmi_path.relative_to(report._coot_root)),
        "mmdb_exists":  mmdb_path.exists(),
        "gemmi_exists": gemmi_path.exists(),
    }


@router.post("/api/tests/write")
def write_test_files(req: WriteTestRequest):
    """Write MMDB and/or Gemmi test .cc files to disk."""
    written: dict[str, str] = {}
    for variant, content in [("mmdb", req.mmdb_test), ("gemmi", req.gemmi_test)]:
        if not content.strip():
            continue
        path = get_test_file_path(req.rel_source_path, req.fn_name, variant)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(wrap_test_content(content, variant, req.fn_name))
        written[variant] = str(path.relative_to(report._coot_root))
    return {"ok": True, "written": written}


@router.post("/api/tests/compile-run")
async def compile_and_run(req: CompileRunRequest):
    """Write the test file, compile it, run it, stream all output as plain text."""
    test_file = get_test_file_path(req.rel_source_path, req.fn_name, req.variant)
    test_file.parent.mkdir(parents=True, exist_ok=True)
    test_file.write_text(wrap_test_content(req.test_code, req.variant, req.fn_name))

    rel_path   = test_file.relative_to(report._coot_root)
    bin_dir    = test_file.parent / "bin"
    bin_dir.mkdir(exist_ok=True)
    output_bin = bin_dir / test_file.stem

    compile_cmd = make_compile_cmd(test_file, output_bin)
    if args.cxx_flags:
        compile_cmd += f" {args.cxx_flags}"

    async def stream_output():
        yield f"[WRITE]   {rel_path}\n"
        yield f"[COMPILE] {compile_cmd}\n\n"
        yield f"[WORKING_DIR] {str(report._coot_root)}\n\n"

        proc = await asyncio.create_subprocess_shell(
            compile_cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            cwd=str(report._coot_root),
        )
        assert proc.stdout is not None
        async for raw in proc.stdout:
            yield raw.decode(errors="replace")
        await proc.wait()

        if proc.returncode != 0:
            yield f"\n[COMPILE FAILED]  exit code {proc.returncode}\n"
            return

        yield f"\n[COMPILE OK]\n[RUN]     {output_bin.relative_to(report._coot_root)}\n\n"

        run_proc = await asyncio.create_subprocess_exec(
            str(output_bin),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
        )
        assert run_proc.stdout is not None
        async for raw in run_proc.stdout:
            yield raw.decode(errors="replace")
        await run_proc.wait()

        code = run_proc.returncode
        yield f"\n{'[PASSED]' if code == 0 else f'[FAILED]  exit code {code}'}\n"

    return StreamingResponse(stream_output(), media_type="text/plain")


# ── Git commit ────────────────────────────────────────────────────────────────

@router.post("/api/tests/git-commit")
async def git_commit_test(req: GitCommitRequest):
    """Git-add and commit test files."""
    variants = ["mmdb", "gemmi"] if req.variant == "both" else [req.variant]
    files_to_add = [
        str(get_test_file_path(req.rel_source_path, req.fn_name, v))
        for v in variants
        if get_test_file_path(req.rel_source_path, req.fn_name, v).exists()
    ]

    async def stream_output():
        if not files_to_add:
            yield "[GIT ERROR] No existing test files found to commit — write them first\n"
            return

        rel_paths = [str(Path(f).relative_to(report._coot_root)) for f in files_to_add]
        yield f"[GIT ADD]  {', '.join(rel_paths)}\n"

        add_proc = await asyncio.create_subprocess_shell(
            "git add " + " ".join(f'"{f}"' for f in files_to_add),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            cwd=str(report._coot_root),
        )
        assert add_proc.stdout is not None
        async for raw in add_proc.stdout:
            yield raw.decode(errors="replace")
        await add_proc.wait()

        if add_proc.returncode != 0:
            yield f"[GIT ADD FAILED]  exit code {add_proc.returncode}\n"
            return

        safe_msg = req.commit_message.replace('"', "'")
        commit_proc = await asyncio.create_subprocess_shell(
            f'git commit -m "{safe_msg}"',
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            cwd=str(report._coot_root),
        )
        assert commit_proc.stdout is not None
        async for raw in commit_proc.stdout:
            yield raw.decode(errors="replace")
        await commit_proc.wait()

        if commit_proc.returncode != 0:
            yield f"\n[GIT COMMIT FAILED]  exit code {commit_proc.returncode}\n"
        else:
            yield "\n[COMMITTED]\n"

    return StreamingResponse(stream_output(), media_type="text/plain")


# ── Validate & fix (single test, single variant) ─────────────────────────────

@router.post("/api/tests/validate-fix")
async def validate_and_fix(req: ValidateFixRequest):
    """Compile + run a test variant; on failure, ask the LLM to fix it.
    Streams the events emitted by `compile_run_fix_loop`, mapped to the
    NDJSON shape the existing frontend expects.
    """
    async def stream():
        path = get_test_file_path(req.rel_source_path, req.fn_name, req.variant)

        yield _ndjson({
            "type": "start", "max": MAX_TEST_RETRIES, "variant": req.variant,
        })

        async for evt in compile_run_fix_loop(
            model=req.model,
            fn_name=req.fn_name,
            variant=req.variant,
            test_code=req.test_code,
            test_path=path,
        ):
            t = evt["type"]
            if t == "loop_start":
                continue  # already emitted "start"
            if t == "compiling":
                yield _ndjson({"type": "compiling", "attempt": evt["attempt"], "max": evt["max"]})
            elif t == "compile_output":
                yield _ndjson({"type": "compile_output", "text": evt["text"], "ok": evt["ok"]})
            elif t == "running":
                yield _ndjson({"type": "running", "attempt": evt["attempt"]})
            elif t == "run_output":
                yield _ndjson({"type": "run_output", "text": evt["text"], "ok": evt["ok"]})
            elif t == "fixing":
                yield _ndjson({"type": "fixing", "attempt": evt["attempt"]})
            elif t == "fixed_code":
                yield _ndjson({"type": "fixed_code", "code": evt["code"]})
            elif t == "loop_done":
                yield _ndjson({
                    "type": "done", "status": evt["status"],
                    "attempts": evt["attempts"], "code": evt["code"],
                })

    return StreamingResponse(stream(), media_type="application/x-ndjson")


# ── Batch generation ──────────────────────────────────────────────────────────

def _function_source(rel_path: str, fn_line: int, fn_end: int) -> str:
    """Slice a function's source out of its containing file."""
    full_path = os.path.join(report._coot_root, rel_path)
    with open(full_path, "r", errors="replace") as f:
        lines = f.readlines()
    start = max(0, fn_line - 1)
    end   = min(len(lines), fn_end)
    return "".join(lines[start:end])


@router.post("/api/tests/generate-all")
async def generate_all_tests(req: GenerateAllRequest):
    """Iterate every function, generate both tests, then compile/fix each.

    The work for one function is:
        1. probe the real MMDB function for ground-truth output (oracle)
        2. ask the LLM for both MMDB and Gemmi tests in one shot
        3. parse the response into two variants
        4. for each variant: compile, run, and ask the LLM to fix on failure

    Streams NDJSON: start | progress | skip | oracle_event | oracle |
                    attempt | done | error | finish
    """
    async def stream_progress():
        tests = load_tests()
        work = [
            (entry["rel_path"], fn)
            for entry in report._files
            for fn in entry.get("functions", [])
        ]
        total = len(work)

        yield _ndjson({"type": "start", "total": total})

        done = skipped = errors = 0

        for rel_path, fn in work:
            fn_name: str = fn["name"]
            fn_line: int = fn["line"]
            fn_end:  int = fn["end_line"]
            mmdb_syms    = fn.get("mmdb_symbols", [])
            key = f"{rel_path}::{fn_name}:{fn_line}"

            # ── Skip if both variants are already on disk ────────────────
            if req.skip_existing:
                rec = tests.get(key, {})
                if rec.get("mmdb_test") and rec.get("gemmi_test"):
                    skipped += 1
                    yield _ndjson({
                        "type": "skip", "key": key, "fn": fn_name, "file": rel_path,
                        "done": done, "skipped": skipped, "errors": errors, "total": total,
                    })
                    continue

            yield _ndjson({
                "type": "progress", "key": key, "fn": fn_name, "file": rel_path,
                "done": done, "skipped": skipped, "errors": errors, "total": total,
            })

            try:
                source = _function_source(rel_path, fn_line, fn_end)

                # ── 1. Oracle probe ──────────────────────────────────────
                probe = None
                async for evt in run_oracle_for_function(
                    model=req.model,
                    function_name=fn_name,
                    source_code=source,
                    mmdb_symbols=mmdb_syms,
                    additional_instructions=req.additional_instructions,
                    rel_source_path=rel_path,
                ):
                    if evt["type"] == "pipeline_result":
                        probe = evt["result"]
                        continue
                    # Tag every event with which function it belongs to.
                    yield _ndjson({**evt, "key": key, "fn": fn_name, "file": rel_path})

                yield _ndjson({
                    "type": "oracle", "key": key, "fn": fn_name, "file": rel_path,
                    "ok":       bool(probe and probe.ok),
                    "stage":    probe.stage    if probe else "none",
                    "attempts": probe.attempts if probe else 0,
                    "lines":    (probe.probe_lines[:50] if probe else []),
                })

                oracle_out = oracle_output_text(probe)
                probe_src  = oracle_probe_source(probe)

                # ── 2. Ask the LLM for both tests in one go ──────────────
                prompt = build_generate_test(
                    fn_name, source, mmdb_syms, "both",
                    req.additional_instructions,
                    oracle_output=oracle_out,
                    rel_source_path=rel_path,
                    probe_source=probe_src,
                    probe_pdb_path=PROBE_PDB_PATH,
                )
                raw = await call_ollama(req.model, prompt)
                mmdb_raw, gemmi_raw = parse_both_sections(raw)
                initial = {
                    "mmdb":  strip_fences(mmdb_raw),
                    "gemmi": strip_fences(gemmi_raw),
                }

                # ── 3. Compile/run/fix each variant ──────────────────────
                variant_results: dict[str, dict] = {}
                for variant in ("mmdb", "gemmi"):
                    code = initial[variant]
                    if not code.strip():
                        variant_results[variant] = {
                            "code": "", "status": "skip", "attempts": 0,
                        }
                        continue

                    path = get_test_file_path(rel_path, fn_name, variant)
                    final = {"status": "fail", "attempts": 0, "code": code}

                    async for evt in compile_run_fix_loop(
                        model=req.model,
                        fn_name=fn_name,
                        variant=variant,
                        test_code=code,
                        test_path=path,
                    ):
                        if evt["type"] == "fixing":
                            yield _ndjson({
                                "type": "attempt",
                                "key": key, "fn": fn_name, "file": rel_path,
                                "variant": variant,
                                "attempt": evt["attempt"], "max": MAX_TEST_RETRIES,
                                "error": (evt.get("error") or "")[:600],
                                "done": done, "skipped": skipped, "errors": errors, "total": total,
                            })
                        elif evt["type"] == "loop_done":
                            final = {
                                "status":   evt["status"],
                                "attempts": evt["attempts"],
                                "code":     evt["code"],
                            }

                    variant_results[variant] = final

                # ── 4. Persist whatever we ended up with ─────────────────
                rec = tests.get(key, {"notes": "", "status": "draft"})
                for variant, result in variant_results.items():
                    if result["code"].strip():
                        rec[f"{variant}_test"] = result["code"]
                tests[key] = rec
                save_tests(tests)

                mmdb_r  = variant_results.get("mmdb",  {"status": "skip", "attempts": 0})
                gemmi_r = variant_results.get("gemmi", {"status": "skip", "attempts": 0})
                done += 1
                yield _ndjson({
                    "type": "done", "key": key, "fn": fn_name, "file": rel_path,
                    "done": done, "skipped": skipped, "errors": errors, "total": total,
                    "mmdb_status":   mmdb_r["status"],  "mmdb_attempts":  mmdb_r["attempts"],
                    "gemmi_status":  gemmi_r["status"], "gemmi_attempts": gemmi_r["attempts"],
                })

            except Exception as ex:
                errors += 1
                yield _ndjson({
                    "type": "error", "key": key, "fn": fn_name, "file": rel_path,
                    "done": done, "skipped": skipped, "errors": errors, "total": total,
                    "message": str(ex),
                })

        yield _ndjson({
            "type": "finish", "done": done, "skipped": skipped, "errors": errors, "total": total,
        })

    return StreamingResponse(stream_progress(), media_type="application/x-ndjson")
