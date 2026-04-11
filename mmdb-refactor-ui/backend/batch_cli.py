#!/usr/bin/env python3
"""
Command-line batch test generator.

Runs the same oracle → generate → compile/fix pipeline as the web UI but
prints progress to stdout, so you can leave it running in a terminal.

Usage:
    cd mmdb-refactor-ui/backend
    python batch_cli.py [options]

Examples:
    # Generate MMDB tests, skip functions already done:
    python batch_cli.py --target mmdb

    # Use a specific model, regenerate everything:
    python batch_cli.py --target mmdb --model gemma4:31b --no-skip

    # Generate both variants:
    python batch_cli.py --target both --model gemma4:31b
"""
import argparse
import asyncio
import os
from pathlib import Path

# ── CLI args (parsed before importing config so --report etc. work) ──────────
_cli = argparse.ArgumentParser(
    description="Batch MMDB/Gemmi test generator (CLI)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
_cli.add_argument("--target",  default="mmdb", choices=["mmdb", "gemmi", "both"],
                  help="Which test variant(s) to generate")
_cli.add_argument("--model",   default="gemma4:31b",
                  help="Ollama model name")
_cli.add_argument("--no-skip", action="store_true",
                  help="Regenerate tests even if they already exist")
_cli.add_argument("--instructions", default="",
                  help="Additional instructions passed to every prompt")
_cli_args, _ = _cli.parse_known_args()

# ── Backend imports (config.py parses its own --report / --probe-pdb etc.) ───
import report
import prompts
from config import PROBE_PDB_PATH, PROBE_WORKDIR, REPORT_PATH, args
from ollama import call_ollama
from pipeline import (
    compile_run_fix_loop,
    oracle_output_text,
    oracle_probe_source,
    run_oracle_for_function,
)
from prompts import (
    build_generate_test,
    parse_both_sections,
    strip_fences,
    system_context_for_test_target,
)
from storage import load_tests, save_tests
from test_utils import get_test_file_path


# ── Helpers ───────────────────────────────────────────────────────────────────

def _function_source(rel_path: str, fn_line: int, fn_end: int) -> str:
    full_path = os.path.join(report._coot_root, rel_path)
    with open(full_path, "r", errors="replace") as f:
        lines = f.readlines()
    return "".join(lines[max(0, fn_line - 1):min(len(lines), fn_end)])


def _status_char(status: str) -> str:
    return {"pass": "✓", "fail": "✗", "skip": "·"}.get(status, "?")


def _indent(text: str, prefix: str = "           ") -> str:
    """Indent every line of text for nested output."""
    return "\n".join(prefix + ln for ln in text.splitlines()) if text.strip() else ""


def _print_oracle_event(evt: dict) -> None:
    """Pretty-print a single oracle pipeline event."""
    t = evt.get("type")

    if t == "oracle_event":
        inner = evt.get("event", {})
        inner_t = inner.get("type", "")

        if inner_t == "info":
            print(f"             {inner['text']}")

        elif inner_t == "llm_call":
            print(
                f"             LLM call  attempt={inner['attempt']}/{inner['max']}"
                f"  model={inner['model']}  prompt={inner['prompt_chars']:,} chars"
            )

        elif inner_t == "compile_cmd":
            print(f"             compile:  {inner['cmd']}")

        elif inner_t == "compile_out":
            ok_str = "ok" if inner["ok"] else "FAILED"
            print(f"             compile → {ok_str}  ({inner['elapsed_s']}s)")
            if not inner["ok"] and inner.get("text"):
                print(_indent(inner["text"].strip(), "               "))

        elif inner_t == "run_cmd":
            print(f"             run:      {inner['cmd']}")

        elif inner_t == "run_out":
            ok_str = "ok" if inner["ok"] else "FAILED"
            lines = inner.get("probe_lines", [])
            print(f"             run   → {ok_str}  ({inner['elapsed_s']}s)  {len(lines)} PROBE line(s)")
            for ln in lines[:10]:
                print(f"               {ln}")
            if len(lines) > 10:
                print(f"               … {len(lines) - 10} more")

        elif inner_t == "retry":
            print(f"             retry #{inner['attempt']}  reason: {inner['reason']}")

        elif inner_t == "llm_response":
            print(f"             LLM response  {inner['source_chars']:,} chars  ({inner['elapsed_s']}s)")

    elif t == "oracle_done":
        ok_str = "ok" if evt["ok"] else "FAILED"
        print(
            f"           oracle done → {ok_str}  stage={evt['stage']}"
            f"  attempts={evt['attempts']}"
        )
        prompt_file = PROBE_WORKDIR / "probe_prompt.txt"
        probe_file  = PROBE_WORKDIR / "probe.cc"
        print(f"           probe prompt  → {prompt_file}")
        print(f"           probe source  → {probe_file}")


def _print_fix_event(evt: dict, variant: str) -> None:
    """Pretty-print a single compile/fix loop event."""
    t = evt.get("type")

    if t == "compiling":
        print(f"           [{variant}] compile  attempt {evt['attempt']}/{evt['max']} …", flush=True)

    elif t == "compile_output":
        ok_str = "ok" if evt["ok"] else "FAILED"
        print(f"           [{variant}] compile → {ok_str}")
        if not evt["ok"] and evt.get("text"):
            print(_indent(evt["text"].strip(), "             "))

    elif t == "running":
        print(f"           [{variant}] run …", flush=True)

    elif t == "run_output":
        ok_str = "passed" if evt["ok"] else "FAILED"
        print(f"           [{variant}] run    → {ok_str}")
        if not evt["ok"] and evt.get("text"):
            print(_indent(evt["text"].strip(), "             "))

    elif t == "fixing":
        print(f"           [{variant}] asking LLM to fix (attempt {evt['attempt']}) …", flush=True)


# ── Main batch loop ───────────────────────────────────────────────────────────

async def run_batch(
    target: str,
    model: str,
    skip_existing: bool,
    additional_instructions: str,
) -> None:
    variants_to_run = ("mmdb", "gemmi") if target == "both" else (target,)
    tests = load_tests()

    work = [
        (entry["rel_path"], fn)
        for entry in report._files
        for fn in entry.get("functions", [])
    ]
    total = len(work)
    done = skipped = errors = 0

    print(f"[BATCH] target={target}  model={model}  total={total}  skip_existing={skip_existing}")
    print(f"[BATCH] report  → {REPORT_PATH}")
    print(f"[BATCH] pdb     → {PROBE_PDB_PATH}")
    print(f"[BATCH] workdir → {PROBE_WORKDIR}")
    print()

    for idx, (rel_path, fn) in enumerate(work, 1):
        fn_name: str = fn["name"]
        fn_line: int = fn["line"]
        fn_end:  int = fn["end_line"]
        mmdb_syms    = fn.get("mmdb_symbols", [])
        key = f"{rel_path}::{fn_name}:{fn_line}"
        short = fn_name.split("::")[-1]
        prefix = f"[{idx:>4}/{total}]"

        # ── Skip ──────────────────────────────────────────────────────────
        if skip_existing:
            rec = tests.get(key, {})
            if all(rec.get(f"{v}_test") for v in variants_to_run):
                skipped += 1
                print(f"{prefix} · {short}  (skip)")
                continue

        print(f"{prefix} {fn_name}")
        print(f"         source: {rel_path}:{fn_line}")

        try:
            source = _function_source(rel_path, fn_line, fn_end)

            # ── Oracle ────────────────────────────────────────────────────
            probe = None
            if target in ("mmdb", "both"):
                print(f"         oracle:")
                async for evt in run_oracle_for_function(
                    model=model,
                    function_name=fn_name,
                    source_code=source,
                    mmdb_symbols=mmdb_syms,
                    additional_instructions=additional_instructions,
                    rel_source_path=rel_path,
                ):
                    if evt["type"] == "pipeline_result":
                        probe = evt["result"]
                    else:
                        _print_oracle_event(evt)

            oracle_out = oracle_output_text(probe)
            probe_src  = oracle_probe_source(probe)

            # ── Generate ──────────────────────────────────────────────────
            print(f"         generate ({target}) …", flush=True)
            if target == "both":
                prompt = build_generate_test(
                    fn_name, source, mmdb_syms, "both",
                    additional_instructions,
                    oracle_output=oracle_out,
                    rel_source_path=rel_path,
                    probe_source=probe_src,
                    probe_pdb_path=PROBE_PDB_PATH,
                )
                raw = await call_ollama(model, prompt, system=system_context_for_test_target("both"))
                mmdb_raw, gemmi_raw = parse_both_sections(raw)
                initial = {"mmdb": strip_fences(mmdb_raw), "gemmi": strip_fences(gemmi_raw)}
            else:
                prompt = build_generate_test(
                    fn_name, source, mmdb_syms, target,
                    additional_instructions,
                    oracle_output=oracle_out,
                    rel_source_path=rel_path,
                    probe_source=probe_src,
                    probe_pdb_path=PROBE_PDB_PATH,
                )
                raw = await call_ollama(model, prompt, system=system_context_for_test_target(target))
                initial = {target: strip_fences(raw)}

            # Save generate prompt next to probe prompt for inspection
            gen_prompt_file = PROBE_WORKDIR / f"generate_prompt_{target}.txt"
            gen_prompt_file.write_text(prompt)
            print(f"         generate prompt → {gen_prompt_file}")
            print(f"         generate done  ({sum(len(c) for c in initial.values()):,} chars)")

            # ── Compile / fix ─────────────────────────────────────────────
            variant_results: dict[str, dict] = {}
            for variant in variants_to_run:
                code = initial.get(variant, "")
                if not code.strip():
                    print(f"           [{variant}] empty — skipping")
                    variant_results[variant] = {"code": "", "status": "skip", "attempts": 0}
                    continue

                test_path = get_test_file_path(rel_path, fn_name, variant)
                print(f"           [{variant}] test file → {test_path}")
                final: dict = {"status": "fail", "attempts": 0, "code": code}

                async for evt in compile_run_fix_loop(
                    model=model,
                    fn_name=fn_name,
                    variant=variant,
                    test_code=code,
                    test_path=test_path,
                ):
                    _print_fix_event(evt, variant)
                    if evt["type"] == "loop_done":
                        final = {"status": evt["status"], "attempts": evt["attempts"], "code": evt["code"]}

                ch = _status_char(final["status"])
                print(
                    f"           [{variant}] {ch} {final['status']}"
                    f"  ({final['attempts']} attempt{'s' if final['attempts'] != 1 else ''})"
                )
                variant_results[variant] = final

            # ── Persist ───────────────────────────────────────────────────
            rec = tests.get(key, {"notes": "", "status": "draft"})
            for variant, result in variant_results.items():
                if result["code"].strip():
                    rec[f"{variant}_test"] = result["code"]
            tests[key] = rec
            save_tests(tests)
            done += 1
            print()

        except KeyboardInterrupt:
            print("\n[BATCH] interrupted by user")
            break
        except Exception as ex:
            errors += 1
            print(f"         ERROR: {ex}")
            print()

    print(f"[BATCH] finished — done={done}  skipped={skipped}  errors={errors}  total={total}")


if __name__ == "__main__":
    report.load_and_update(REPORT_PATH)
    prompts.init_docs(args.mmdb_docs)
    asyncio.run(run_batch(
        target=_cli_args.target,
        model=_cli_args.model,
        skip_existing=not _cli_args.no_skip,
        additional_instructions=_cli_args.instructions,
    ))
