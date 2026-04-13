#!/usr/bin/env python3
"""
Summary of generated test coverage.

Reads tests.json and the usage report to show how many functions have
passing/written/missing tests for each variant.

Usage:
    cd mmdb-refactor-ui/backend
    python summary_cli.py [options]

Examples:
    python summary_cli.py
    python summary_cli.py --target mmdb
    python summary_cli.py --verbose
"""
import argparse
import sys

_cli = argparse.ArgumentParser(
    description="Summarise generated test coverage from tests.json",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
_cli.add_argument("--target", default="both", choices=["mmdb", "gemmi", "both"],
                  help="Which variant(s) to report on")
_cli.add_argument("--verbose", "-v", action="store_true",
                  help="List every function and its status")
_cli.add_argument("--missing", action="store_true",
                  help="Only show functions with no test written")
_cli_args, _ = _cli.parse_known_args()

import report
from config import REPORT_PATH, TESTS_PATH, args
from storage import load_tests


def _bar(n: int, total: int, width: int = 30) -> str:
    filled = int(width * n / total) if total else 0
    return "[" + "█" * filled + "·" * (width - filled) + "]"


def main() -> None:
    report.load_and_update(REPORT_PATH)
    tests = load_tests()

    target = _cli_args.target
    variants = ["mmdb", "gemmi"] if target == "both" else [target]

    # Build the full function list from the report
    all_functions: list[tuple[str, str]] = [
        (entry["rel_path"], fn)
        for entry in report._files
        for fn in entry.get("functions", [])
    ]
    total = len(all_functions)

    # Per-variant counters
    counts: dict[str, dict[str, int]] = {
        v: {"written": 0, "not_written": 0} for v in variants
    }
    details: list[tuple[str, str, dict[str, str]]] = []  # (key, fn_name, {variant: status})

    for rel_path, fn in all_functions:
        fn_name: str = fn["name"]
        fn_line: int = fn["line"]
        key = f"{rel_path}::{fn_name}:{fn_line}"
        rec = tests.get(key, {})

        statuses: dict[str, str] = {}
        for v in variants:
            code = rec.get(f"{v}_test", "").strip()
            if code:
                statuses[v] = "written"
                counts[v]["written"] += 1
            else:
                statuses[v] = "not_written"
                counts[v]["not_written"] += 1

        details.append((key, fn_name, statuses))

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\n{'─'*60}")
    print(f"  Test Coverage Summary")
    print(f"  Report : {REPORT_PATH}")
    print(f"  Store  : {TESTS_PATH}")
    print(f"  Total functions: {total}")
    print(f"{'─'*60}")

    for v in variants:
        c = counts[v]
        written = c["written"]
        not_written = c["not_written"]
        pct = 100 * written / total if total else 0
        bar = _bar(written, total)
        print(f"\n  [{v.upper()}]")
        print(f"    Written     : {written:>5}  {bar}  {pct:.1f}%")
        print(f"    Not written : {not_written:>5}")

    print(f"\n{'─'*60}\n")

    # ── Verbose / missing listing ─────────────────────────────────────────────
    if _cli_args.verbose or _cli_args.missing:
        show_missing_only = _cli_args.missing

        for key, fn_name, statuses in details:
            if show_missing_only:
                # Only show if at least one variant is not written
                if all(s == "written" for s in statuses.values()):
                    continue

            tags = "  ".join(
                f"{v.upper()}:{'✓' if s == 'written' else '✗'}"
                for v, s in statuses.items()
            )
            short = fn_name.split("::")[-1]
            print(f"  {tags}  {fn_name}")


if __name__ == "__main__":
    main()
