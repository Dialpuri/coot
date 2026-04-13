#!/usr/bin/env python3
"""
Summary of generated test coverage.

Reads tests.json and the usage report to show how many functions have
passing/written/missing tests for each variant.  With --run it also
compiles and executes every written test to report pass/fail counts.

Usage:
    cd mmdb-refactor-ui/backend
    python summary_cli.py [options]

Examples:
    python summary_cli.py
    python summary_cli.py --target mmdb
    python summary_cli.py --run
    python summary_cli.py --run --verbose
    python summary_cli.py --missing
"""
import argparse
import asyncio
import sys

_cli = argparse.ArgumentParser(
    description="Summarise generated test coverage from tests.json",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
_cli.add_argument("--target", default="both", choices=["mmdb", "gemmi", "both"],
                  help="Which variant(s) to report on")
_cli.add_argument("--run", action="store_true",
                  help="Compile and run every written test to check pass/fail")
_cli.add_argument("--verbose", "-v", action="store_true",
                  help="List every function and its status")
_cli.add_argument("--missing", action="store_true",
                  help="Only show functions with no test written")
_cli_args, _ = _cli.parse_known_args()

import report
from config import REPORT_PATH, TESTS_PATH, args
from compiler import compile_test, run_test
from storage import load_tests
from test_utils import get_test_file_path


def _bar(n: int, total: int, width: int = 30) -> str:
    filled = int(width * n / total) if total else 0
    return "[" + "█" * filled + "·" * (width - filled) + "]"


# Status values used in details list
# "not_written" | "pass" | "compile_error" | "fail"

_STATUS_CHAR = {
    "not_written":   "·",
    "pass":          "✓",
    "compile_error": "C",
    "fail":          "✗",
}


async def _run_existing(rel_path: str, fn_name: str, variant: str) -> str:
    """Run an already-compiled test binary; return status string."""
    test_path = get_test_file_path(rel_path, fn_name, variant)
    bin_path = test_path.parent / "bin" / test_path.stem
    if not bin_path.exists():
        return "no_binary"
    ok, _ = await run_test(test_path)
    return "pass" if ok else "fail"


async def main() -> None:
    report.load_and_update(REPORT_PATH)
    tests = load_tests()

    target = _cli_args.target
    variants = ["mmdb", "gemmi"] if target == "both" else [target]
    do_run = _cli_args.run

    # Build the full function list from the report
    all_functions: list[tuple[str, str]] = [
        (entry["rel_path"], fn)
        for entry in report._files
        for fn in entry.get("functions", [])
    ]
    total = len(all_functions)

    # Per-variant counters
    counts: dict[str, dict[str, int]] = {
        v: {"not_written": 0, "written": 0, "pass": 0, "no_binary": 0, "fail": 0}
        for v in variants
    }
    details: list[tuple[str, str, str, dict[str, str]]] = []  # (key, fn_name, rel_path, {v: status})

    # Gather compile/run tasks if --run
    if do_run:
        tasks: dict[tuple[str, str], asyncio.Task] = {}
        for entry in report._files:
            rel_path = entry["rel_path"]
            for fn in entry.get("functions", []):
                fn_name = fn["name"]
                fn_line = fn["line"]
                key = f"{rel_path}::{fn_name}:{fn_line}"
                rec = tests.get(key, {})
                for v in variants:
                    code = rec.get(f"{v}_test", "").strip()
                    if code:
                        tasks[(key, v)] = asyncio.create_task(
                            _run_existing(rel_path, fn_name, v)
                        )

        written_count = len(tasks) // len(variants) if variants else 0
        print(f"  Running {len(tasks)} test(s) …", flush=True)
        await asyncio.gather(*tasks.values())
        run_results: dict[tuple[str, str], str] = {k: t.result() for k, t in tasks.items()}
    else:
        run_results = {}

    for entry in report._files:
        rel_path = entry["rel_path"]
        for fn in entry.get("functions", []):
            fn_name: str = fn["name"]
            fn_line: int = fn["line"]
            key = f"{rel_path}::{fn_name}:{fn_line}"
            rec = tests.get(key, {})

            statuses: dict[str, str] = {}
            for v in variants:
                code = rec.get(f"{v}_test", "").strip()
                if not code:
                    status = "not_written"
                elif do_run:
                    status = run_results.get((key, v), "fail")
                else:
                    status = "written"
                statuses[v] = status
                counts[v][status if do_run else ("written" if code else "not_written")] += 1

            details.append((key, fn_name, rel_path, statuses))

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\n{'─'*60}")
    print(f"  Test Coverage Summary")
    print(f"  Report : {REPORT_PATH}")
    print(f"  Store  : {TESTS_PATH}")
    print(f"  Total functions: {total}")
    print(f"{'─'*60}")

    for v in variants:
        c = counts[v]
        print(f"\n  [{v.upper()}]")
        if do_run:
            passed      = c["pass"]
            no_binary   = c["no_binary"]
            failed      = c["fail"]
            not_written = c["not_written"]
            written     = passed + no_binary + failed
            pct_written = 100 * written / total if total else 0
            pct_pass    = 100 * passed   / total if total else 0
            print(f"    Written     : {written:>5}  {_bar(written, total)}  {pct_written:.1f}%")
            print(f"    Passing     : {passed:>5}  {_bar(passed, total)}  {pct_pass:.1f}%")
            print(f"    Failing     : {failed:>5}")
            print(f"    No binary   : {no_binary:>5}")
            print(f"    Not written : {not_written:>5}")
        else:
            written     = c["written"]
            not_written = c["not_written"]
            pct = 100 * written / total if total else 0
            print(f"    Written     : {written:>5}  {_bar(written, total)}  {pct:.1f}%")
            print(f"    Not written : {not_written:>5}")

    print(f"\n{'─'*60}\n")

    # ── Verbose / missing listing ─────────────────────────────────────────────
    if _cli_args.verbose or _cli_args.missing:
        show_missing_only = _cli_args.missing

        for key, fn_name, rel_path, statuses in details:
            if show_missing_only:
                if all(s not in ("not_written",) for s in statuses.values()):
                    continue

            tags = "  ".join(
                f"{v.upper()}:{_STATUS_CHAR.get(s, '?')}"
                for v, s in statuses.items()
            )
            print(f"  {tags}  {fn_name}")


if __name__ == "__main__":
    asyncio.run(main())
