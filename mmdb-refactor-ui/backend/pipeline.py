"""
Shared building blocks for the test-generation endpoints.

The HTTP endpoints in `routers/tests.py` are mostly thin wrappers around
three pipelines that all need to compose:

    1. run_oracle_for_function   — drive the LLM probe loop and capture
                                   ground-truth output values from a real run.
    2. stream_test_generation    — stream a Google Test from Ollama, optionally
                                   seeded with the oracle's measured values.
    3. compile_run_fix_loop      — write a test to disk, compile it, run it,
                                   and ask the LLM to fix it on failure.

Splitting these out keeps `routers/tests.py` focused on HTTP concerns and
makes the data flow explicit:

    oracle  →  prompt with measured values  →  generated test  →  fix loop

Every helper here is an async generator that yields event dicts. The events
are forwarded by the HTTP layer as NDJSON lines so the frontend can show
progress in real time.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import AsyncIterator

from compiler import compile_test, run_test
from config import MAX_TEST_RETRIES
from ollama import call_ollama, stream_ollama
from oracle import ProbeResult, run_mmdb_oracle
from prompts import (
    PROBE_SYSTEM_CONTEXT,
    STRATEGY_CONTEXT,
    build_fix,
    build_generate_test,
    strip_fences,
)
from test_utils import wrap_test_content


# ── 1. Oracle ─────────────────────────────────────────────────────────────────

async def run_oracle_for_function(
    *,
    model: str,
    function_name: str,
    source_code: str,
    mmdb_symbols: list[str],
    additional_instructions: str,
    rel_source_path: str,
) -> AsyncIterator[dict]:
    """Wrap `run_mmdb_oracle` so it always yields one final `oracle_done`
    event followed by a `pipeline_result` carrying the `ProbeResult`.

    Yields:
        oracle_start                 — emitted before the probe loop
        oracle_event {event: ...}    — every event from `run_mmdb_oracle`
        oracle_done  {ok, stage, attempts, probe_lines, raw, probe_source}
        pipeline_result {result: ProbeResult | None}
    """
    yield {"type": "oracle_start"}

    result: ProbeResult | None = None
    try:
        async for evt in run_mmdb_oracle(
            model=model,
            function_name=function_name,
            source_code=source_code,
            mmdb_symbols=mmdb_symbols,
            additional_instructions=additional_instructions,
            rel_source_path=rel_source_path,
        ):
            if evt.get("type") == "done":
                result = evt["result"]
                continue
            yield {"type": "oracle_event", "event": evt}
    except Exception as ex:
        yield {"type": "oracle_skip", "reason": f"oracle crashed: {ex}"}

    if result is not None:
        yield {
            "type":         "oracle_done",
            "ok":           result.ok,
            "stage":        result.stage,
            "attempts":     result.attempts,
            "probe_lines":  result.probe_lines,
            "raw":          result.raw_stdout[:4000],
            "probe_source": result.probe_source,
        }

    yield {"type": "pipeline_result", "result": result}


def oracle_output_text(result: ProbeResult | None) -> str:
    """Convert a successful ProbeResult into the multi-line string the
    test-generation prompt expects, or "" if the probe failed."""
    if result is None or not result.ok:
        return ""
    return "\n".join(result.probe_lines)


def oracle_probe_source(result: ProbeResult | None) -> str:
    """Return the working probe `main()` source from a successful run, or "".

    Passed into the test-generation prompt so the LLM has a verbatim,
    already-compiled example of how to load the data and call the function
    instead of inventing its own (broken) loader.
    """
    if result is None or not result.ok:
        return ""
    return result.probe_source or ""


# ── 2. Test generation ────────────────────────────────────────────────────────

async def stream_test_generation(
    *,
    model: str,
    function_name: str,
    source_code: str,
    mmdb_symbols: list[str],
    target: str,
    additional_instructions: str,
    oracle_output: str,
    rel_source_path: str,
    probe_source: str = "",
    probe_pdb_path: str = "",
) -> AsyncIterator[dict]:
    """Stream a test from Ollama, surfacing thinking + response chunks
    separately so the UI can show the model's reasoning.

    `probe_source` / `probe_pdb_path` come from a successful oracle run and
    are injected into the prompt so the LLM uses the working loader as the
    basis for the test instead of inventing one.

    Yields:
        test_start                            — first event
        test_thinking {text}                  — reasoning chunks (thinking models only)
        test_chunk    {text}                  — response chunks (the actual code)
        test_done     {full_text}             — emitted once on success
        error         {message}               — only on failure
    """
    prompt = build_generate_test(
        function_name, source_code, mmdb_symbols,
        target, additional_instructions,
        oracle_output=oracle_output,
        rel_source_path=rel_source_path,
        probe_source=probe_source,
        probe_pdb_path=probe_pdb_path,
    )

    yield {"type": "test_start"}

    parts: list[str] = []
    try:
        async for evt in stream_ollama(model, prompt, system=STRATEGY_CONTEXT):
            if evt["kind"] == "thinking":
                yield {"type": "test_thinking", "text": evt["text"]}
            else:
                parts.append(evt["text"])
                yield {"type": "test_chunk", "text": evt["text"]}
    except Exception as ex:
        yield {"type": "error", "message": str(ex)}
        return

    yield {"type": "test_done", "full_text": "".join(parts)}


# ── 3. Compile / run / fix loop ───────────────────────────────────────────────

async def compile_run_fix_loop(
    *,
    model: str,
    fn_name: str,
    variant: str,                # "mmdb" | "gemmi"
    test_code: str,
    test_path: Path,
    max_retries: int = MAX_TEST_RETRIES,
) -> AsyncIterator[dict]:
    """Iterate compile → (run on success) → ask LLM to fix on failure.

    Yields:
        loop_start      {variant, max}
        compiling       {attempt, max}
        compile_output  {ok, text}
        running         {attempt}
        run_output      {ok, text}
        fixing          {attempt, error}        — error is the captured failure text
        fixed_code      {code}
        loop_done       {status: "pass"|"fail", attempts, code}
    """
    test_path.parent.mkdir(parents=True, exist_ok=True)
    current_code = test_code

    yield {"type": "loop_start", "variant": variant, "max": max_retries}

    for attempt in range(1, max_retries + 1):
        test_path.write_text(wrap_test_content(current_code, variant, fn_name))

        yield {"type": "compiling", "attempt": attempt, "max": max_retries}
        compile_ok, compile_out = await compile_test(test_path)
        yield {"type": "compile_output", "ok": compile_ok, "text": compile_out}

        if compile_ok:
            yield {"type": "running", "attempt": attempt}
            run_ok, run_out = await run_test(test_path)
            yield {"type": "run_output", "ok": run_ok, "text": run_out}
            if run_ok:
                yield {
                    "type": "loop_done", "status": "pass",
                    "attempts": attempt, "code": current_code,
                }
                return
            error_out = f"Compile: OK\nRun failed:\n{run_out}"
        else:
            error_out = compile_out

        if attempt >= max_retries:
            yield {
                "type": "loop_done", "status": "fail",
                "attempts": attempt, "code": current_code,
            }
            return

        # Ask the LLM to repair the test using the captured error.
        yield {"type": "fixing", "attempt": attempt, "error": error_out}
        fix_prompt = build_fix(fn_name, variant, current_code, error_out, attempt)
        fixed_raw  = await call_ollama(model, fix_prompt)
        current_code = strip_fences(fixed_raw) or current_code
        yield {"type": "fixed_code", "code": current_code}
