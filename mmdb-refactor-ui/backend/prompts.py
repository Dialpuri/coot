"""
Prompt loading and construction.

Templates live in the prompts/ directory as plain-text files.
They use Python's string.Template syntax: $variable or ${variable}.
C++ curly braces in the template body are safe — no escaping needed.
"""
import re
from pathlib import Path
from string import Template

import doc_extractor
from config import MAX_TEST_RETRIES

_PROMPTS_DIR = Path(__file__).parent / "prompts"


def _load(name: str) -> Template:
    return Template((_PROMPTS_DIR / name).read_text())


# ── Static context (loaded once) ─────────────────────────────────────────────
STRATEGY_CONTEXT: str = (_PROMPTS_DIR / "strategy_context.txt").read_text()

# ── Templates (loaded once) ───────────────────────────────────────────────────
_T_REFACTOR          = _load("refactor.txt")
_T_TEST_MMDB         = _load("generate_test_mmdb.txt")
_T_TEST_GEMMI        = _load("generate_test_gemmi.txt")
_T_TEST_BOTH         = _load("generate_test_both.txt")
_T_FIX               = _load("fix_test.txt")

# ── MMDB docs (injected at startup via init_docs) ─────────────────────────────
_mmdb_docs_markdown: str = ""


def init_docs(path: str) -> None:
    """Load MMDB API docs once at server startup."""
    global _mmdb_docs_markdown
    _mmdb_docs_markdown = doc_extractor.load_docs(path)
    if _mmdb_docs_markdown:
        print(f"Loaded MMDB docs from {path} ({len(_mmdb_docs_markdown):,} chars)")
    else:
        print(f"WARNING: MMDB docs not found at {path} — prompts will lack API context")


def _api_context(symbols: list[str]) -> str:
    """Return a markdown snippet with only the MMDB doc sections relevant to symbols."""
    if not _mmdb_docs_markdown or not symbols:
        return ""
    excerpt = doc_extractor.extract_for_symbols(_mmdb_docs_markdown, symbols)
    if not excerpt:
        return ""
    return (
        "\n\n## Relevant MMDB API Reference\n\n"
        "The following is extracted from the MMDB API documentation. "
        "Use it to understand the types and methods in the source code below.\n\n"
        + excerpt
    )


def _additional(text: str) -> str:
    """Format optional additional instructions line, or empty string."""
    s = text.strip()
    return f"Additional instructions: {s}\n" if s else ""


# ── Public prompt builders ────────────────────────────────────────────────────

def build_refactor(function_name: str, code: str, mmdb_symbols: list[str],
                   additional_instructions: str = "") -> str:
    return _T_REFACTOR.substitute(
        function_name=function_name,
        symbols_list=", ".join(mmdb_symbols) if mmdb_symbols else "none listed",
        additional_instructions=_additional(additional_instructions),
        code=code,
    )


def build_generate_test(function_name: str, source_code: str, mmdb_symbols: list[str],
                        target: str, additional_instructions: str = "") -> str:
    """Build a test-generation prompt for target='mmdb'|'gemmi'|'both'."""
    common = dict(
        function_name=function_name,
        symbols_list=", ".join(mmdb_symbols) if mmdb_symbols else "none",
        additional_instructions=_additional(additional_instructions),
        api_context=_api_context(mmdb_symbols),
        source_code=source_code,
    )
    if target == "mmdb":
        return _T_TEST_MMDB.substitute(**common)
    if target == "gemmi":
        return _T_TEST_GEMMI.substitute(**common)
    return _T_TEST_BOTH.substitute(**common)


def build_fix(fn_name: str, variant: str, current_code: str,
              error: str, attempt: int) -> str:
    if variant == "mmdb":
        inc_hint = "Required headers:\n#include <gtest/gtest.h>\n#include <mmdb2/mmdb_manager.h>"
    else:
        inc_hint = ("Required headers:\n#include <gtest/gtest.h>\n"
                    "#include <gemmi/structure.hpp>\n#include <gemmi/model.hpp>")
    return _T_FIX.substitute(
        fn_name=fn_name,
        variant=variant.upper(),
        attempt=attempt,
        max_retries=MAX_TEST_RETRIES,
        inc_hint=inc_hint,
        current_code=current_code,
        error=error[:2000],
    )


# ── Parse helpers (used by batch generation) ─────────────────────────────────

def parse_both_sections(raw: str) -> tuple[str, str]:
    mmdb_marker  = "=== MMDB TEST ==="
    gemmi_marker = "=== GEMMI TEST ==="
    mi, gi = raw.find(mmdb_marker), raw.find(gemmi_marker)
    if mi == -1 and gi == -1:
        return raw, ""
    if mi != -1 and gi == -1:
        return raw[mi + len(mmdb_marker):].strip(), ""
    if mi == -1:
        return "", raw[gi + len(gemmi_marker):].strip()
    return raw[mi + len(mmdb_marker):gi].strip(), raw[gi + len(gemmi_marker):].strip()


def strip_fences(code: str) -> str:
    m = re.match(r'^```(?:cpp)?\s*\n([\s\S]*?)```\s*$', code.strip())
    return m.group(1) if m else code
