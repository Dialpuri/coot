"""Endpoints: /api/tests/*, /api/generate-test"""
import json
import os

import httpx
from fastapi import APIRouter, Query
from fastapi.responses import StreamingResponse
from pathlib import Path

import report
from compiler import compile_test, run_test, make_compile_cmd
from config import OLLAMA_URL, MAX_TEST_RETRIES, args
from models import (
    TestRecord, GenerateTestRequest, WriteTestRequest,
    CompileRunRequest, GitCommitRequest, GenerateAllRequest, ValidateFixRequest,
)
from ollama import call_ollama
from prompts import (
    STRATEGY_CONTEXT, build_generate_test, build_fix,
    parse_both_sections, strip_fences,
)
from storage import load_tests, save_tests
from test_utils import get_test_file_path, wrap_test_content

router = APIRouter()


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
    """Return the test record for a specific function, falling back to .cc files on disk."""
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


# ── Generate (streaming) ──────────────────────────────────────────────────────

@router.post("/api/generate-test")
async def generate_test(req: GenerateTestRequest):
    prompt = build_generate_test(
        req.function_name, req.source_code, req.mmdb_symbols,
        req.target, req.additional_instructions,
    )
    payload = {
        "model":  req.model,
        "system": STRATEGY_CONTEXT,
        "prompt": prompt,
        "stream": True,
    }

    async def stream_ollama():
        async with httpx.AsyncClient(timeout=300.0) as client:
            async with client.stream("POST", OLLAMA_URL, json=payload) as response:
                if response.status_code != 200:
                    body = await response.aread()
                    yield f"Error from Ollama ({response.status_code}): {body.decode()}"
                    return
                async for line in response.aiter_lines():
                    if line.strip():
                        try:
                            data = json.loads(line)
                            chunk = data.get("response", "")
                            if chunk:
                                yield chunk
                        except json.JSONDecodeError:
                            pass

    return StreamingResponse(stream_ollama(), media_type="text/plain")


# ── File paths / write / compile-run ─────────────────────────────────────────

@router.get("/api/tests/file-paths")
def get_test_file_paths(rel_source_path: str = Query(...), fn_name: str = Query(...)):
    """Return where MMDB and Gemmi test files would be written."""
    mmdb_path  = get_test_file_path(rel_source_path, fn_name, "mmdb")
    gemmi_path = get_test_file_path(rel_source_path, fn_name, "gemmi")
    return {
        "mmdb":        str(mmdb_path.relative_to(report._coot_root)),
        "gemmi":       str(gemmi_path.relative_to(report._coot_root)),
        "mmdb_exists": mmdb_path.exists(),
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
    """Write the test file, compile it, run it, stream all output."""
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
        import asyncio
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
        import asyncio
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


# ── Validate & fix ────────────────────────────────────────────────────────────

@router.post("/api/tests/validate-fix")
async def validate_and_fix(req: ValidateFixRequest):
    """
    Compile + run a single test variant. On failure, ask the LLM to fix it and
    retry up to MAX_TEST_RETRIES times. Streams NDJSON progress events.

    Event types: start | compiling | compile_output | running | run_output |
                 fixing | fixed_code | done
    """
    async def stream():
        path = get_test_file_path(req.rel_source_path, req.fn_name, req.variant)
        path.parent.mkdir(parents=True, exist_ok=True)
        current_code = req.test_code

        yield json.dumps({"type": "start", "max": MAX_TEST_RETRIES, "variant": req.variant}) + "\n"

        for attempt in range(1, MAX_TEST_RETRIES + 1):
            path.write_text(wrap_test_content(current_code, req.variant, req.fn_name))

            yield json.dumps({"type": "compiling", "attempt": attempt, "max": MAX_TEST_RETRIES}) + "\n"
            compile_ok, compile_out = await compile_test(path)
            yield json.dumps({"type": "compile_output", "text": compile_out, "ok": compile_ok}) + "\n"

            if compile_ok:
                yield json.dumps({"type": "running", "attempt": attempt}) + "\n"
                run_ok, run_out = await run_test(path)
                yield json.dumps({"type": "run_output", "text": run_out, "ok": run_ok}) + "\n"
                if run_ok:
                    yield json.dumps({"type": "done", "status": "pass", "attempts": attempt, "code": current_code}) + "\n"
                    return
                error_out = f"Compile: OK\nRun failed:\n{run_out}"
            else:
                error_out = compile_out

            if attempt < MAX_TEST_RETRIES:
                yield json.dumps({"type": "fixing", "attempt": attempt}) + "\n"
                fix_prompt  = build_fix(req.fn_name, req.variant, current_code, error_out, attempt)
                fixed_raw   = await call_ollama(req.model, fix_prompt)
                current_code = strip_fences(fixed_raw) or current_code
                yield json.dumps({"type": "fixed_code", "code": current_code}) + "\n"
            else:
                yield json.dumps({"type": "done", "status": "fail", "attempts": attempt, "code": current_code}) + "\n"

    return StreamingResponse(stream(), media_type="application/x-ndjson")


# ── Batch generation ──────────────────────────────────────────────────────────

@router.post("/api/tests/generate-all")
async def generate_all_tests(req: GenerateAllRequest):
    """
    Iterate every function, generate both tests, compile+run each, retry up to
    MAX_TEST_RETRIES times on failure using the compile error as LLM feedback.
    Streams NDJSON: start | progress | skip | attempt | done | error | finish
    """
    async def stream_progress():
        tests = load_tests()
        work = [
            (entry["rel_path"], fn)
            for entry in report._files
            for fn in entry.get("functions", [])
        ]
        total = len(work)
        yield json.dumps({"type": "start", "total": total}) + "\n"

        done = skipped = errors = 0

        for rel_path, fn in work:
            fn_name: str = fn["name"]
            fn_line: int = fn["line"]
            fn_end: int  = fn["end_line"]
            key = f"{rel_path}::{fn_name}:{fn_line}"

            if req.skip_existing:
                rec = tests.get(key, {})
                if rec.get("mmdb_test") and rec.get("gemmi_test"):
                    skipped += 1
                    yield json.dumps({
                        "type": "skip", "key": key, "fn": fn_name, "file": rel_path,
                        "done": done, "skipped": skipped, "errors": errors, "total": total,
                    }) + "\n"
                    continue

            yield json.dumps({
                "type": "progress", "key": key, "fn": fn_name, "file": rel_path,
                "done": done, "skipped": skipped, "errors": errors, "total": total,
            }) + "\n"

            try:
                full_path = os.path.join(report._coot_root, rel_path)
                with open(full_path, "r", errors="replace") as f:
                    all_lines = f.readlines()
                source = "".join(all_lines[max(0, fn_line - 1):min(len(all_lines), fn_end)])

                prompt  = build_generate_test(fn_name, source, fn.get("mmdb_symbols", []),
                                              "both", req.additional_instructions)
                raw     = await call_ollama(req.model, prompt)
                mmdb_raw, gemmi_raw = parse_both_sections(raw)
                initial = {"mmdb": strip_fences(mmdb_raw), "gemmi": strip_fences(gemmi_raw)}

                variant_results: dict[str, dict] = {}
                for variant in ("mmdb", "gemmi"):
                    current_code = initial[variant]
                    if not current_code.strip():
                        variant_results[variant] = {"code": "", "status": "skip", "attempts": 0}
                        continue

                    path = get_test_file_path(rel_path, fn_name, variant)
                    path.parent.mkdir(parents=True, exist_ok=True)

                    for attempt in range(1, MAX_TEST_RETRIES + 1):
                        path.write_text(wrap_test_content(current_code, variant, fn_name))
                        compile_ok, compile_out = await compile_test(path)

                        if compile_ok:
                            run_ok, run_out = await run_test(path)
                            if run_ok:
                                variant_results[variant] = {"code": current_code, "status": "pass", "attempts": attempt}
                                break
                            error_out = f"Compile: OK\nRun output:\n{run_out}"
                        else:
                            error_out = compile_out

                        if attempt < MAX_TEST_RETRIES:
                            yield json.dumps({
                                "type": "attempt",
                                "key": key, "fn": fn_name, "file": rel_path,
                                "variant": variant, "attempt": attempt, "max": MAX_TEST_RETRIES,
                                "error": error_out[:600],
                                "done": done, "skipped": skipped, "errors": errors, "total": total,
                            }) + "\n"
                            fix_prompt   = build_fix(fn_name, variant, current_code, error_out, attempt)
                            fixed_raw    = await call_ollama(req.model, fix_prompt)
                            current_code = strip_fences(fixed_raw) or current_code
                        else:
                            variant_results[variant] = {"code": current_code, "status": "fail", "attempts": attempt}

                rec = tests.get(key, {"notes": "", "status": "draft"})
                for variant, result in variant_results.items():
                    if result["code"].strip():
                        rec[f"{variant}_test"] = result["code"]
                tests[key] = rec
                save_tests(tests)

                mmdb_r  = variant_results.get("mmdb",  {"status": "skip", "attempts": 0})
                gemmi_r = variant_results.get("gemmi", {"status": "skip", "attempts": 0})
                done += 1
                yield json.dumps({
                    "type": "done", "key": key, "fn": fn_name, "file": rel_path,
                    "done": done, "skipped": skipped, "errors": errors, "total": total,
                    "mmdb_status":  mmdb_r["status"],  "mmdb_attempts":  mmdb_r["attempts"],
                    "gemmi_status": gemmi_r["status"], "gemmi_attempts": gemmi_r["attempts"],
                }) + "\n"

            except Exception as ex:
                errors += 1
                yield json.dumps({
                    "type": "error", "key": key, "fn": fn_name, "file": rel_path,
                    "done": done, "skipped": skipped, "errors": errors, "total": total,
                    "message": str(ex),
                }) + "\n"

        yield json.dumps({
            "type": "finish", "done": done, "skipped": skipped, "errors": errors, "total": total,
        }) + "\n"

    return StreamingResponse(stream_progress(), media_type="application/x-ndjson")
