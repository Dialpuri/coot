"""
Oracle probes: ask the LLM for a standalone C++ main() that runs a function on
a real PDB file, compile it, run it, and capture its stdout so the measured
output values can be fed back into the test-generation prompt.

`run_mmdb_oracle` is an async generator that yields progress events describing
each step (LLM call, compile command, run output, retries, etc.). The final
event is always `{"type": "done", "result": ProbeResult}`.
"""
import asyncio
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import AsyncIterator

import report
from compiler import make_probe_compile_cmd
from config import PROBE_PDB_PATH, MAX_PROBE_RETRIES, DEV_STREAM_LLM, args
from ollama import call_ollama, stream_ollama
from prompts import PROBE_SYSTEM_CONTEXT, build_probe_mmdb, strip_fences


@dataclass
class ProbeResult:
    ok: bool
    raw_stdout: str              # everything the probe printed (or compile error)
    probe_lines: list[str]       # just the "PROBE: ..." lines
    probe_source: str            # the C++ source that was compiled
    attempts: int
    stage: str                   # "llm" | "compile" | "run" | "ok"


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


async def run_mmdb_oracle(
    model: str,
    function_name: str,
    source_code: str,
    mmdb_symbols: list[str],
    additional_instructions: str = "",
    pdb_path: str | None = None,
    rel_source_path: str = "",
) -> AsyncIterator[dict]:
    """
    Async generator. Yields event dicts describing each step of the probe
    pipeline. The last event is always:

        {"type": "done", "result": ProbeResult}

    Event types (all include "type"):
      info          — a human-readable status line  {"text": str}
      llm_call      — about to call Ollama          {"attempt": int, "max": int, "model": str}
      llm_response  — LLM returned probe source     {"attempt": int, "source": str, "elapsed_s": float}
      compile_cmd   — compile command about to run  {"cmd": str}
      compile_out   — compiler stdout/stderr        {"ok": bool, "text": str, "elapsed_s": float}
      run_cmd       — probe binary about to run     {"cmd": str, "pdb": str}
      run_out       — probe stdout                  {"ok": bool, "text": str, "probe_lines": list, "elapsed_s": float}
      retry         — about to retry                {"attempt": int, "reason": str}
      done          — final result                  {"result": ProbeResult}
    """
    pdb = pdb_path or PROBE_PDB_PATH
    yield {"type": "info", "text": f"probe PDB: {pdb}"}
    if not Path(pdb).exists():
        yield {"type": "done", "result": ProbeResult(
            ok=False,
            raw_stdout=f"probe PDB file not found: {pdb}",
            probe_lines=[], probe_source="", attempts=0, stage="run",
        )}
        return

    base_prompt = build_probe_mmdb(
        function_name, source_code, mmdb_symbols, additional_instructions,
        rel_source_path=rel_source_path,
    )
    prompt = base_prompt
    last_source = ""
    last_error = ""
    yield {"type": "info", "text": f"base_prompt → {base_prompt}"}

    with tempfile.TemporaryDirectory(prefix="mmdb_probe_") as tmp:
        tmp_dir = Path(tmp)
        probe_src = tmp_dir / "probe.cc"
        probe_bin = tmp_dir / "probe"
        yield {"type": "info", "text": f"probe workdir: {tmp_dir}"}

        for attempt in range(1, MAX_PROBE_RETRIES + 1):
            # ── 1. Ask the LLM for probe source ──────────────────────────────
            yield {
                "type": "llm_call",
                "attempt": attempt,
                "max": MAX_PROBE_RETRIES,
                "model": model,
                "prompt_chars": len(prompt),
                "streaming": DEV_STREAM_LLM,
            }
            t0 = time.monotonic()

            if DEV_STREAM_LLM:
                # Dev mode: stream chunks live into the terminal so we can
                # watch the LLM write the probe. Same total compute as the
                # non-streaming path; only the transport differs.
                raw_parts: list[str] = []
                try:
                    async for chunk in stream_ollama(model, prompt, system=PROBE_SYSTEM_CONTEXT):
                        raw_parts.append(chunk)
                        yield {"type": "llm_chunk", "attempt": attempt, "text": chunk}
                except Exception as ex:
                    yield {"type": "info", "text": f"LLM stream error: {ex}"}
                raw = "".join(raw_parts)
            else:
                raw = await call_ollama(model, prompt, system=PROBE_SYSTEM_CONTEXT)

            llm_elapsed = time.monotonic() - t0
            source = strip_fences(raw).strip()

            yield {
                "type": "llm_response",
                "attempt": attempt,
                "source": source,
                "source_chars": len(source),
                "elapsed_s": round(llm_elapsed, 2),
            }

            if not source:
                last_error = "LLM returned empty probe source"
                yield {"type": "info", "text": "LLM returned empty source — retrying"}
                prompt = f"{base_prompt}\n\nPrevious attempt produced no code. Try again."
                if attempt < MAX_PROBE_RETRIES:
                    yield {"type": "retry", "attempt": attempt, "reason": "empty LLM response"}
                continue

            # Pre-compile sanity check: the probe must be pure MMDB. If the LLM
            # slipped gemmi in anyway (ignoring both system and user prompts),
            # reject and retry with explicit feedback instead of wasting a compile.
            if "gemmi" in source.lower():
                last_source = source
                last_error = "probe source contained gemmi — rejected before compile"
                yield {
                    "type": "info",
                    "text": "probe contained 'gemmi' — rejecting and retrying (MMDB only)",
                }
                if attempt < MAX_PROBE_RETRIES:
                    yield {"type": "retry", "attempt": attempt, "reason": "contained gemmi"}
                    prompt = (
                        f"{base_prompt}\n\n"
                        "Your previous probe contained `gemmi` code. That is forbidden. "
                        "This probe MUST be pure MMDB2 — no `gemmi::` types, no "
                        "`#include <gemmi/...>`, no Gemmi helpers anywhere. "
                        "Rewrite the probe using only MMDB2 and the real coot header "
                        "for the function under test.\n\n"
                        f"Previous (rejected) source:\n```cpp\n{source}\n```"
                    )
                continue  # on the last attempt this exits the loop → final done

            last_source = source
            probe_src.write_text(source)
            yield {"type": "info", "text": f"prompt → {prompt}"}

            # ── 2. Compile ───────────────────────────────────────────────────
            cmd = _make_compile_cmd(probe_src, probe_bin)
            yield {"type": "compile_cmd", "cmd": cmd}
            t0 = time.monotonic()
            compile_ok, compile_out = await _compile_probe(cmd)
            compile_elapsed = time.monotonic() - t0
            yield {
                "type": "compile_out",
                "ok": compile_ok,
                "text": compile_out,
                "elapsed_s": round(compile_elapsed, 2),
            }

            if not compile_ok:
                last_error = compile_out
                if attempt < MAX_PROBE_RETRIES:
                    yield {"type": "retry", "attempt": attempt, "reason": "compile failed"}
                    prompt = (
                        f"{base_prompt}\n\n"
                        f"Your previous probe failed to compile. Error:\n```\n{compile_out[:2000]}\n```\n"
                        f"Previous source:\n```cpp\n{source}\n```\n"
                        "Fix the problem and output the full corrected probe source."
                    )
                    continue
                yield {"type": "done", "result": ProbeResult(
                    ok=False, raw_stdout=compile_out, probe_lines=[],
                    probe_source=source, attempts=attempt, stage="compile",
                )}
                return

            # ── 3. Run ───────────────────────────────────────────────────────
            yield {"type": "run_cmd", "cmd": f"{probe_bin} {pdb}", "pdb": pdb}
            t0 = time.monotonic()
            run_ok, run_out = await _run_probe(probe_bin, pdb)
            run_elapsed = time.monotonic() - t0
            probe_lines = _extract_probe_lines(run_out)
            yield {
                "type": "run_out",
                "ok": run_ok,
                "text": run_out,
                "probe_lines": probe_lines,
                "elapsed_s": round(run_elapsed, 2),
            }

            if run_ok and probe_lines:
                yield {"type": "done", "result": ProbeResult(
                    ok=True, raw_stdout=run_out, probe_lines=probe_lines,
                    probe_source=source, attempts=attempt, stage="ok",
                )}
                return

            last_error = run_out
            if attempt < MAX_PROBE_RETRIES:
                reason = "no PROBE: lines emitted" if run_ok else "probe exited non-zero"
                yield {"type": "retry", "attempt": attempt, "reason": reason}
                hint = ("No PROBE: lines were printed — remember every output line must start with 'PROBE: '."
                        if run_ok else
                        f"The probe exited with an error. Output:\n```\n{run_out[:2000]}\n```")
                prompt = (
                    f"{base_prompt}\n\n"
                    f"Your previous probe compiled but did not produce usable output. {hint}\n"
                    f"Previous source:\n```cpp\n{source}\n```\n"
                    "Fix the problem and output the full corrected probe source."
                )

        yield {"type": "done", "result": ProbeResult(
            ok=False, raw_stdout=last_error, probe_lines=[],
            probe_source=last_source, attempts=MAX_PROBE_RETRIES, stage="run",
        )}
