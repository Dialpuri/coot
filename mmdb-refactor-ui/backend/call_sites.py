"""
call_sites.py — find real call sites for MMDB symbols in the coot tree.

Why this exists
---------------
The LLM keeps calling non-static MMDB methods statically — `mmdb::Manager::GetModel(1)`
instead of `mgr->GetModel(1)`. The doc extractor tells the model what methods
exist but never shows it the *receiver*. This module fills that gap by greppping
the coot source tree for verbatim call sites and injecting them into the prompt.

There is no LLM in this pipeline. The whole point is that everything we
return came from a real translation unit in coot, so the model can copy it.

Public API
----------
    find_call_sites(symbols, max_per_symbol=3) -> dict[symbol, list[CallSite]]
    format_call_sites_for_prompt(symbols)      -> str   (markdown section)
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

import report
from doc_extractor import _parse_symbol


# ── Tunables ─────────────────────────────────────────────────────────────────

_SOURCE_EXTS = (".cc", ".cpp", ".cxx", ".hh", ".h", ".hpp", ".hxx")

# Directories we never want example code from. mmdb-recon and mmdb-refactor-ui
# are *this* tooling tree (their examples would be circular). build dirs are
# generated. tests/ would surface our own previously-generated tests, which
# defeats the point of grounding in pre-existing real usage.
_EXCLUDE_DIR_NAMES = {
    ".git", "node_modules", ".venv", "venv", "__pycache__",
    "mmdb-recon", "mmdb-refactor-ui", "macos-build",
    "tests",
}
_EXCLUDE_DIR_PREFIXES = ("build",)

# Lines longer than this are usually generated tables, long string literals,
# or noise — not the kind of "see, this is how you call it" example we want.
_MAX_LINE_LENGTH = 220
_DEFAULT_MAX_PER_SYMBOL = 3


# ── Data type ────────────────────────────────────────────────────────────────

@dataclass
class CallSite:
    rel_path: str
    lineno:   int
    text:     str        # the source line, stripped


# ── Module-level caches ──────────────────────────────────────────────────────

_file_list_cache: list[Path] | None = None
_call_site_cache: dict[str, list[CallSite]] = {}


def _gather_source_files() -> list[Path]:
    """Walk the coot tree once and remember every C++ source / header path."""
    global _file_list_cache
    if _file_list_cache is not None:
        return _file_list_cache
    root = report._coot_root
    if not root:
        _file_list_cache = []
        return _file_list_cache
    files: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        # prune walk in-place
        dirnames[:] = [
            d for d in dirnames
            if d not in _EXCLUDE_DIR_NAMES
            and not d.startswith(".")
            and not any(d.startswith(p) for p in _EXCLUDE_DIR_PREFIXES)
        ]
        for fn in filenames:
            if fn.endswith(_SOURCE_EXTS):
                files.append(Path(dirpath) / fn)
    _file_list_cache = files
    print(f"call_sites: indexed {len(files)} source files under {root}")
    return files


# ── Pattern building ─────────────────────────────────────────────────────────

def _pattern_for_symbol(raw: str) -> re.Pattern | None:
    """Build a regex that matches a real call site for the given mmdb symbol.

    For methods (`mmdb::Chain::GetChainID`) we look for `chain->GetChainID(`
    or `chain.GetChainID(` style calls — those show the receiver, which is
    exactly the information the model is missing.

    For bare types (`mmdb::Chain`) we look for variable declarations and
    parameter spellings — `mmdb::Chain *chain` or `mmdb::Chain&`.
    """
    class_name, method_name = _parse_symbol(raw)
    if not class_name and not method_name:
        return None
    if method_name:
        # Receiver-based call: ->Method(  or .Method(
        # Word boundary on the left so we don't match SubMethod inside BigMethod.
        return re.compile(rf"(?:->|\.)\s*{re.escape(method_name)}\s*\(")
    # Type usage: mmdb::Class followed by space/star/amp/end-of-token
    return re.compile(rf"\bmmdb::{re.escape(class_name)}\b")


# ── Result selection ─────────────────────────────────────────────────────────

def _normalise(text: str) -> str:
    """Collapse whitespace so two visually-identical hits dedupe."""
    return re.sub(r"\s+", " ", text).strip()


def _select_diverse(hits: list[CallSite], max_n: int) -> list[CallSite]:
    """Pick up to `max_n` distinct, illustrative hits.

    Strategy:
      1. Drop dupes by normalised text.
      2. Prefer shorter lines (clearer examples).
      3. Prefer hits from .cc files over headers (real call sites, not declarations).
    """
    seen_text: set[str] = set()
    deduped: list[CallSite] = []
    for h in hits:
        key = _normalise(h.text)
        if key in seen_text:
            continue
        seen_text.add(key)
        deduped.append(h)

    def score(h: CallSite) -> tuple[int, int]:
        is_header = 0 if h.rel_path.endswith((".cc", ".cpp", ".cxx")) else 1
        return (is_header, len(h.text))

    deduped.sort(key=score)
    return deduped[:max_n]


# ── Public API ───────────────────────────────────────────────────────────────

def find_call_sites(
    symbols: list[str],
    max_per_symbol: int = _DEFAULT_MAX_PER_SYMBOL,
) -> dict[str, list[CallSite]]:
    """Return up to `max_per_symbol` real call sites per symbol.

    Symbols already in `_call_site_cache` are served from memory. Cache misses
    trigger one combined pass over the source tree, which sets the cache for
    every requested symbol so a future call with overlapping symbols is free.
    """
    results: dict[str, list[CallSite]] = {}
    pending: dict[str, re.Pattern] = {}

    for raw in symbols:
        if raw in _call_site_cache:
            results[raw] = _call_site_cache[raw]
            continue
        pat = _pattern_for_symbol(raw)
        if pat is None:
            _call_site_cache[raw] = []
            results[raw] = []
            continue
        pending[raw] = pat

    if not pending:
        return results

    files = _gather_source_files()
    if not files:
        for sym in pending:
            _call_site_cache[sym] = []
            results[sym] = []
        return results

    root = Path(report._coot_root)
    per_symbol_hits: dict[str, list[CallSite]] = {s: [] for s in pending}

    # Cap how many raw hits we collect per symbol before _select_diverse runs;
    # otherwise a method like `GetChainID` would scan the whole tree even
    # though we only need a handful of distinct examples.
    raw_cap = max_per_symbol * 8
    saturated: set[str] = set()

    for path in files:
        if len(saturated) == len(pending):
            break
        try:
            with path.open("r", errors="replace") as fh:
                lines = fh.readlines()
        except OSError:
            continue
        try:
            rel = str(path.relative_to(root))
        except ValueError:
            rel = str(path)

        for lineno, line in enumerate(lines, 1):
            if len(line) > _MAX_LINE_LENGTH:
                continue
            stripped = line.strip()
            if not stripped or stripped.startswith("//") or stripped.startswith("*"):
                continue
            for sym, pat in pending.items():
                if sym in saturated:
                    continue
                if pat.search(line):
                    bucket = per_symbol_hits[sym]
                    bucket.append(CallSite(rel_path=rel, lineno=lineno, text=stripped))
                    if len(bucket) >= raw_cap:
                        saturated.add(sym)

    for sym, hits in per_symbol_hits.items():
        chosen = _select_diverse(hits, max_per_symbol)
        _call_site_cache[sym] = chosen
        results[sym] = chosen

    return results


def format_call_sites_for_prompt(
    symbols: list[str],
    max_per_symbol: int = _DEFAULT_MAX_PER_SYMBOL,
) -> str:
    """Render real call sites as a prompt section, or '' if none were found.

    The output deliberately includes file:line annotations so the LLM can
    see these are real places in the codebase, not synthesised examples.
    """
    if not symbols:
        return ""

    hits = find_call_sites(symbols, max_per_symbol=max_per_symbol)
    blocks: list[str] = []
    for sym in symbols:
        sym_hits = hits.get(sym, [])
        if not sym_hits:
            continue
        lines = [f"// {h.rel_path}:{h.lineno}\n{h.text}" for h in sym_hits]
        blocks.append(f"### `{sym}`\n```cpp\n" + "\n\n".join(lines) + "\n```")

    if not blocks:
        return ""

    body = "\n\n".join(blocks)
    return (
        "\n\n## Real call sites in the coot source tree\n\n"
        "Below are verbatim lines from the existing coot codebase showing how "
        "each MMDB symbol is actually used — including the **receiver** "
        "(`chain->GetChainID()`, not `mmdb::Chain::GetChainID()`). Mirror "
        "these patterns in your output. If a method is called on an instance "
        "in every example, it is an instance method — do NOT call it via the "
        "class name.\n\n"
        f"{body}\n"
    )
