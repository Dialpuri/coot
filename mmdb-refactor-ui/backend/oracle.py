"""
Oracle probe pipeline.

Goal
----
Test data must be grounded in real execution, not in the LLM's guesses about
what a function returns. The oracle does that:

    1. Ask the LLM to write a tiny standalone C++ probe — a `main()` that
       loads a real PDB file, calls the unmodified MMDB function, and prints
       its output prefixed with "PROBE: ".
    2. Compile the probe.
    3. Run it on the configured PDB.
    4. Capture the "PROBE: ..." lines from stdout.
    5. Hand those lines back to the test-generation prompt as ground-truth
       expected values.

If any step fails (LLM produced gemmi code, didn't compile, ran but emitted
no PROBE lines, etc.) we feed the error back to the LLM and retry, up to
`MAX_PROBE_RETRIES` times.

Public API
----------
`run_mmdb_oracle(...)` is an async generator that yields progress events for
the UI. The very last event is always:

    {"type": "done", "result": ProbeResult}

Event vocabulary
----------------
    info          {"text": str}                                 — status line
    llm_call      {"attempt": int, "max": int, "model": str,
                   "prompt_chars": int, "streaming": bool}      — about to call Ollama
    llm_thinking  {"attempt": int, "text": str}                 — reasoning chunk (streaming)
    llm_chunk     {"attempt": int, "text": str}                 — response chunk (streaming)
    llm_response  {"attempt": int, "source": str,
                   "source_chars": int, "elapsed_s": float}     — full LLM output
    compile_cmd   {"cmd": str}                                  — about to compile
    compile_out   {"ok": bool, "text": str, "elapsed_s": float} — compiler output
    run_cmd       {"cmd": str, "pdb": str}                      — about to run
    run_out       {"ok": bool, "text": str,
                   "probe_lines": list, "elapsed_s": float}     — probe stdout
    retry         {"attempt": int, "reason": str}               — going around again
    done          {"result": ProbeResult}                       — final outcome
"""
import asyncio
import time
from dataclasses import dataclass
from pathlib import Path
from typing import AsyncIterator

import report
from compiler import make_probe_compile_cmd
from config import MAX_PROBE_RETRIES, PROBE_PDB_PATH, PROBE_WORKDIR, args
from ollama import stream_ollama
from prompts import PROBE_SYSTEM_CONTEXT, build_probe_mmdb, strip_fences
from test_utils import sanitize_fn_name


# ── Result type ───────────────────────────────────────────────────────────────

@dataclass
class ProbeResult:
    ok: bool
    raw_stdout: str          # everything the probe printed (or compile error)
    probe_lines: list[str]   # only the "PROBE: ..." lines
    probe_source: str        # the C++ source that was compiled
    attempts: int
    stage: str               # "llm" | "compile" | "run" | "ok"


# ── Subprocess helpers ────────────────────────────────────────────────────────

def _make_compile_cmd(probe_src: Path, output_bin: Path) -> str:
    cmd = make_probe_compile_cmd(probe_src, output_bin)
    if args.cxx_flags:
        cmd += f" {args.cxx_flags}"
    return cmd


async def _compile_probe(cmd: str) -> tuple[bool, str]:
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
        cwd=str(report._coot_root),
    )
    out, _ = await proc.communicate()
    return proc.returncode == 0, out.decode(errors="replace")


async def _run_probe(output_bin: Path, pdb_path: str) -> tuple[bool, str]:
    proc = await asyncio.create_subprocess_exec(
        str(output_bin), pdb_path,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
    )
    out, _ = await proc.communicate()
    return proc.returncode == 0, out.decode(errors="replace")


def _extract_probe_lines(stdout: str) -> list[str]:
    return [ln.strip() for ln in stdout.splitlines() if ln.strip().startswith("PROBE:")]


# ── LLM step ──────────────────────────────────────────────────────────────────

async def _ask_llm_for_probe(
    model: str,
    prompt: str,
    attempt: int,
) -> AsyncIterator[dict]:
    """Stream a probe `main()` from the LLM.

    The probe is the most important LLM call in the whole pipeline — it
    grounds every downstream test in real measured values — so we ALWAYS
    stream it. That keeps the response chunks visible in the UI as they
    arrive AND surfaces any `thinking` chunks from a reasoning model so we
    can watch it work through the API before it commits to code.

    Yields zero or more `llm_thinking` / `llm_chunk` events while the
    model talks, then exactly one `llm_response` event carrying the
    fence-stripped probe source.
    """
    yield {
        "type":         "llm_call",
        "attempt":      attempt,
        "max":          MAX_PROBE_RETRIES,
        "model":        model,
        "prompt_chars": len(prompt),
        "streaming":    True,
    }

    t0 = time.monotonic()
    parts: list[str] = []
    try:
        async for evt in stream_ollama(model, prompt, system=PROBE_SYSTEM_CONTEXT):
            if evt["kind"] == "thinking":
                yield {"type": "llm_thinking", "attempt": attempt, "text": evt["text"]}
            else:
                parts.append(evt["text"])
                yield {"type": "llm_chunk", "attempt": attempt, "text": evt["text"]}
    except Exception as ex:
        yield {"type": "info", "text": f"LLM stream error: {ex}"}

    source = strip_fences("".join(parts)).strip()
    yield {
        "type":         "llm_response",
        "attempt":      attempt,
        "source":       source,
        "source_chars": len(source),
        "elapsed_s":    round(time.monotonic() - t0, 2),
    }


# ── Sanity & retry-prompt builders ────────────────────────────────────────────

def _validate_probe_source(source: str) -> str | None:
    """Return a rejection reason, or None if the source is acceptable."""
    if not source:
        return "empty"
    if "gemmi" in source.lower():
        return "contains_gemmi"
    return None


def _retry_prompt_for_rejection(base_prompt: str, reason: str, source: str) -> str:
    if reason == "empty":
        return f"{base_prompt}\n\nPrevious attempt produced no code. Try again."
    if reason == "contains_gemmi":
        return (
            f"{base_prompt}\n\n"
            "Your previous probe contained `gemmi` code. That is forbidden. "
            "This probe MUST be pure MMDB2 — no `gemmi::` types, no "
            "`#include <gemmi/...>`, no Gemmi helpers anywhere. Rewrite the "
            "probe using only MMDB2 and the real coot header for the function "
            "under test.\n\n"
            f"Previous (rejected) source:\n```cpp\n{source}\n```"
        )
    return base_prompt


def _retry_prompt_for_compile_error(base_prompt: str, compile_out: str, source: str) -> str:
    return (
        f"{base_prompt}\n\n"
        f"Your previous probe failed to compile. Error:\n```\n{compile_out[:2000]}\n```\n"
        f"Previous source:\n```cpp\n{source}\n```\n"
        "Fix the problem and output the full corrected probe source."
    )


def _retry_prompt_for_run_failure(
    base_prompt: str, source: str, run_out: str, run_ok: bool
) -> str:
    if run_ok:
        hint = (
            "No PROBE: lines were printed.\n\n"
            "If the function returns a non-void value, make sure you capture it "
            "and print it as `PROBE: return_value=<value>`.\n\n"
            "If the function returns void, you must instead print the observable "
            "side-effects: read back relevant state from the object(s) that were "
            "passed to the function (e.g. atom counts, coordinates, residue names) "
            "AFTER the call and prefix every printed line with `PROBE: `. "
            "If the function has no observable side-effect, print exactly:\n"
            "    PROBE: void_function=ran_without_error\n"
            "so the oracle knows the call completed."
        )
    else:
        hint = f"The probe exited with an error. Output:\n```\n{run_out[:2000]}\n```"
    return (
        f"{base_prompt}\n\n"
        f"Your previous probe compiled but did not produce usable output. {hint}\n"
        f"Previous source:\n```cpp\n{source}\n```\n"
        "Fix the problem and output the full corrected probe source."
    )


# ── Main pipeline ─────────────────────────────────────────────────────────────

async def run_mmdb_oracle(
    model: str,
    function_name: str,
    source_code: str,
    mmdb_symbols: list[str],
    additional_instructions: str = "",
    pdb_path: str | None = None,
    rel_source_path: str = "",
) -> AsyncIterator[dict]:
    """Drive the LLM → compile → run loop. See module docstring for details."""
    pdb = pdb_path or PROBE_PDB_PATH
    yield {"type": "info", "text": f"probe PDB: {pdb}"}
    if not Path(pdb).exists():
        yield {"type": "done", "result": ProbeResult(
            ok=False,
            raw_stdout=f"probe PDB file not found: {pdb}",
            probe_lines=[], probe_source="", attempts=0, stage="run",
        )}
        return

    # Each function gets its own subdirectory so probes don't overwrite each other.
    fn_workdir = PROBE_WORKDIR / sanitize_fn_name(function_name)
    fn_workdir.mkdir(parents=True, exist_ok=True)
    probe_src = fn_workdir / "probe.cc"
    probe_bin = fn_workdir / "probe"
    yield {"type": "info", "text": f"probe workdir: {fn_workdir}"}

    # If the binary already exists, skip LLM + compile and re-run directly.
    if probe_bin.exists() and probe_src.exists():
        yield {"type": "info", "text": "probe binary found — skipping LLM + compile"}
        yield {"type": "run_cmd", "cmd": f"{probe_bin} {pdb}", "pdb": pdb}
        t0 = time.monotonic()
        run_ok, run_out = await _run_probe(probe_bin, pdb)
        probe_lines = _extract_probe_lines(run_out)
        yield {
            "type":        "run_out",
            "ok":          run_ok,
            "text":        run_out,
            "probe_lines": probe_lines,
            "elapsed_s":   round(time.monotonic() - t0, 2),
        }
        if run_ok and probe_lines:
            yield {"type": "done", "result": ProbeResult(
                ok=True, raw_stdout=run_out, probe_lines=probe_lines,
                probe_source=probe_src.read_text(errors="replace"),
                attempts=0, stage="ok",
            )}
            return
        yield {"type": "info", "text": "cached probe failed — regenerating"}

    base_prompt = build_probe_mmdb(
        function_name, source_code, mmdb_symbols, additional_instructions,
        rel_source_path=rel_source_path, pdb_path=pdb,
    )
    # Persist the prompt next to the probe so devs can re-run it manually.
    (fn_workdir / "probe_prompt.txt").write_text(base_prompt)

    prompt = base_prompt
    last_source = ""
    last_error = ""

    for attempt in range(1, MAX_PROBE_RETRIES + 1):

        # ── 1. Ask the LLM ────────────────────────────────────────────────
        source = ""
        async for evt in _ask_llm_for_probe(model, prompt, attempt):
            if evt["type"] == "llm_response":
                source = evt["source"]
            yield evt

        # ── 2. Sanity-check before wasting a compile ──────────────────────
        rejection = _validate_probe_source(source)
        if rejection is not None:
            last_source = source
            last_error = f"probe rejected before compile: {rejection}"
            yield {"type": "info", "text": f"rejected ({rejection}) — retrying"}
            if attempt < MAX_PROBE_RETRIES:
                yield {"type": "retry", "attempt": attempt, "reason": rejection}
                prompt = _retry_prompt_for_rejection(base_prompt, rejection, source)
            continue

        last_source = source
        probe_src.write_text(source)

        # ── 3. Compile ────────────────────────────────────────────────────
        cmd = _make_compile_cmd(probe_src, probe_bin)
        yield {"type": "compile_cmd", "cmd": cmd}

        t0 = time.monotonic()
        compile_ok, compile_out = await _compile_probe(cmd)
        yield {
            "type":      "compile_out",
            "ok":        compile_ok,
            "text":      compile_out,
            "elapsed_s": round(time.monotonic() - t0, 2),
        }

        if not compile_ok:
            last_error = compile_out
            if attempt < MAX_PROBE_RETRIES:
                yield {"type": "retry", "attempt": attempt, "reason": "compile failed"}
                prompt = _retry_prompt_for_compile_error(base_prompt, compile_out, source)
                continue
            yield {"type": "done", "result": ProbeResult(
                ok=False, raw_stdout=compile_out, probe_lines=[],
                probe_source=source, attempts=attempt, stage="compile",
            )}
            return

        # ── 4. Run ────────────────────────────────────────────────────────
        yield {"type": "run_cmd", "cmd": f"{probe_bin} {pdb}", "pdb": pdb}

        t0 = time.monotonic()
        run_ok, run_out = await _run_probe(probe_bin, pdb)
        probe_lines = _extract_probe_lines(run_out)
        yield {
            "type":        "run_out",
            "ok":          run_ok,
            "text":        run_out,
            "probe_lines": probe_lines,
            "elapsed_s":   round(time.monotonic() - t0, 2),
        }

        if run_ok and probe_lines:
            yield {"type": "done", "result": ProbeResult(
                ok=True, raw_stdout=run_out, probe_lines=probe_lines,
                probe_source=source, attempts=attempt, stage="ok",
            )}
            return

        # ── 5. Run failed → feed error back into prompt and retry ─────────
        last_error = run_out
        if attempt < MAX_PROBE_RETRIES:
            reason = "no PROBE: lines emitted" if run_ok else "probe exited non-zero"
            yield {"type": "retry", "attempt": attempt, "reason": reason}
            prompt = _retry_prompt_for_run_failure(base_prompt, source, run_out, run_ok)

    # Fell off the end of the retry loop without ever succeeding.
    yield {"type": "done", "result": ProbeResult(
        ok=False, raw_stdout=last_error, probe_lines=[],
        probe_source=last_source, attempts=MAX_PROBE_RETRIES, stage="run",
    )}
