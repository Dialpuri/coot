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
import sys
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
from config import PROBE_PDB_PATH, REPORT_PATH, args
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


def _function_source(rel_path: str, fn_line: int, fn_end: int) -> str:
    full_path = os.path.join(report._coot_root, rel_path)
    with open(full_path, "r", errors="replace") as f:
        lines = f.readlines()
    return "".join(lines[max(0, fn_line - 1):min(len(lines), fn_end)])


def _status_char(status: str) -> str:
    return {"pass": "✓", "fail": "✗", "skip": "·"}.get(status, "?")


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

        print(f"{prefix}   {short}", flush=True)

        try:
            source = _function_source(rel_path, fn_line, fn_end)

            # ── Oracle ────────────────────────────────────────────────────
            probe = None
            if target in ("mmdb", "both"):
                print(f"         oracle …", end=" ", flush=True)
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
                ok_str = "ok" if (probe and probe.ok) else "failed"
                attempts = probe.attempts if probe else 0
                print(f"{ok_str} ({attempts} attempt{'s' if attempts != 1 else ''})")

            oracle_out = oracle_output_text(probe)
            probe_src  = oracle_probe_source(probe)

            # ── Generate ──────────────────────────────────────────────────
            print(f"         generate …", end=" ", flush=True)
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
            print("done")

            # ── Compile / fix ─────────────────────────────────────────────
            variant_results: dict[str, dict] = {}
            for variant in variants_to_run:
                code = initial.get(variant, "")
                if not code.strip():
                    variant_results[variant] = {"code": "", "status": "skip", "attempts": 0}
                    continue

                path = get_test_file_path(rel_path, fn_name, variant)
                final: dict = {"status": "fail", "attempts": 0, "code": code}

                print(f"         compile/fix [{variant}] …", end=" ", flush=True)
                async for evt in compile_run_fix_loop(
                    model=model,
                    fn_name=fn_name,
                    variant=variant,
                    test_code=code,
                    test_path=path,
                ):
                    if evt["type"] == "fixing":
                        print(f"retry {evt['attempt']} …", end=" ", flush=True)
                    elif evt["type"] == "loop_done":
                        final = {"status": evt["status"], "attempts": evt["attempts"], "code": evt["code"]}

                ch = _status_char(final["status"])
                print(f"{ch} {final['status']} ({final['attempts']} attempt{'s' if final['attempts'] != 1 else ''})")
                variant_results[variant] = final

            # ── Persist ───────────────────────────────────────────────────
            rec = tests.get(key, {"notes": "", "status": "draft"})
            for variant, result in variant_results.items():
                if result["code"].strip():
                    rec[f"{variant}_test"] = result["code"]
            tests[key] = rec
            save_tests(tests)
            done += 1

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
