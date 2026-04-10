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

# System prompt used for probe generation. The normal STRATEGY_CONTEXT pushes
# the LLM to rewrite everything in Gemmi, which is the opposite of what a probe
# needs — a probe must call the ORIGINAL MMDB function verbatim so it can
# capture ground-truth output. This alternate system prompt keeps the LLM
# focused purely on MMDB.
PROBE_SYSTEM_CONTEXT: str = (
    "You are an expert C++ developer writing a small instrumentation program "
    "that calls an EXISTING MMDB2-based function and prints its real output.\n\n"
    "CRITICAL RULES:\n"
    "- Use MMDB2 types ONLY. Never emit `gemmi::`, `#include <gemmi/...>`, or "
    "any Gemmi code. This is not a refactor — the probe must call the "
    "original unmodified MMDB function.\n"
    "- Do NOT rewrite, port, or reimplement the function. Include its real "
    "header and call it directly.\n"
    "- Do NOT mock, stub, or substitute the function under test.\n"
    "- Load molecular data by having mmdb::Manager read a real PDB file — "
    "never hand-build synthetic mmdb structures.\n"
)

# ── Templates (loaded once) ───────────────────────────────────────────────────
_T_REFACTOR          = _load("refactor.txt")
_T_TEST_MMDB         = _load("generate_test_mmdb.txt")
_T_TEST_GEMMI        = _load("generate_test_gemmi.txt")
_T_TEST_BOTH         = _load("generate_test_both.txt")
_T_FIX               = _load("fix_test.txt")
_T_PROBE_MMDB        = _load("probe_mmdb.txt")

# ── Style guide (injected into test prompts) ──────────────────────────────────
GTEST_STYLE: str = (_PROMPTS_DIR / "google_test.txt").read_text()

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


# ── Header discovery ──────────────────────────────────────────────────────────

_HEADER_EXTS = (".hh", ".h", ".hpp", ".hxx")


def find_header_for_source(rel_source_path: str) -> str:
    """
    Given a .cc path relative to the coot root (e.g. 'coot-utils/coot-coord-utils.cc'),
    return the matching header path (e.g. 'coot-utils/coot-coord-utils.hh') or ''.

    Looks for files with the same stem and a .hh/.h/.hpp/.hxx extension in the
    same directory as the .cc file.
    """
    import report
    if not rel_source_path:
        return ""
    src = Path(report._coot_root) / rel_source_path
    stem_path = src.with_suffix("")
    for ext in _HEADER_EXTS:
        cand = stem_path.with_suffix(ext)
        if cand.exists():
            try:
                return str(cand.relative_to(report._coot_root))
            except ValueError:
                return str(cand)
    return ""


def _format_header_include(rel_source_path: str) -> str:
    """Render a '## Where the function lives' prompt section, or empty string.

    Tells the LLM explicitly which header to `#include` and stresses that the
    function is real — not to be mocked.
    """
    if not rel_source_path:
        return ""
    header = find_header_for_source(rel_source_path)
    if not header:
        return (
            "\n\n## Where the function lives\n\n"
            f"The function is defined in `{rel_source_path}`. Do NOT mock it, "
            "reimplement it, or declare a local stub — `#include` the real "
            "declaration from its header (in the same directory as the .cc "
            "file) and call it directly. The `-I<coot root>` include path is "
            "already on the compile command.\n"
        )
    return (
        "\n\n## Where the function lives\n\n"
        f"The function is declared in `{header}` (defined in `{rel_source_path}`).\n\n"
        f"You MUST include this header and call the real function — do NOT mock it, "
        "reimplement it, or declare a local stub:\n\n"
        f"```cpp\n#include \"{header}\"\n```\n\n"
        "The coot source tree is already on the compile command's include path "
        "(`-I<coot root>`), so this include will resolve without any extra flags.\n"
    )


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
                        target: str, additional_instructions: str = "",
                        oracle_output: str = "",
                        rel_source_path: str = "") -> str:
    """Build a test-generation prompt for target='mmdb'|'gemmi'|'both'.

    If `oracle_output` is provided, it is the captured stdout of a probe run
    that executed the function on a real PDB file. These values are injected
    into the prompt so the LLM can write assertions against real, measured
    outputs instead of guessing.

    If `rel_source_path` is provided, we look up the matching header and tell
    the LLM to `#include` it so it calls the real function instead of mocking.
    """
    common = dict(
        function_name=function_name,
        symbols_list=", ".join(mmdb_symbols) if mmdb_symbols else "none",
        additional_instructions=_additional(additional_instructions),
        api_context=_api_context(mmdb_symbols),
        source_code=source_code,
        gtest_style=GTEST_STYLE,
        oracle_values=_format_oracle(oracle_output),
        header_include=_format_header_include(rel_source_path),
    )
    if target == "mmdb":
        return _T_TEST_MMDB.substitute(**common)
    if target == "gemmi":
        return _T_TEST_GEMMI.substitute(**common)
    return _T_TEST_BOTH.substitute(**common)


def build_probe_mmdb(function_name: str, source_code: str,
                     mmdb_symbols: list[str], additional_instructions: str = "",
                     rel_source_path: str = "") -> str:
    """Build a prompt asking the LLM to write an MMDB probe main()."""
    return _T_PROBE_MMDB.substitute(
        function_name=function_name,
        symbols_list=", ".join(mmdb_symbols) if mmdb_symbols else "none",
        additional_instructions=_additional(additional_instructions),
        api_context=_api_context(mmdb_symbols),
        source_code=source_code,
        header_include=_format_header_include(rel_source_path),
    )


def _format_oracle(raw: str) -> str:
    """Render captured probe output as a prompt section, or empty string."""
    s = raw.strip()
    if not s:
        return ""
    lines = [ln for ln in s.splitlines() if ln.strip().startswith("PROBE:")]
    if not lines:
        return ""
    body = "\n".join(ln.strip() for ln in lines)
    return (
        "\n\n## Measured output from running the function on a real PDB file\n\n"
        "The following values were obtained by compiling and running a probe "
        "program that called this function on a real PDB structure. Use these "
        "as the expected values in your test assertions — do NOT guess or "
        "hardcode different numbers.\n\n"
        f"```\n{body}\n```\n"
    )


def build_fix(fn_name: str, variant: str, current_code: str,
              error: str, attempt: int) -> str:
    if variant == "mmdb":
        inc_hint = "Required headers:\n#include <gtest/gtest.h>\n#include <mmdb2/mmdb_manager.h>\n#include <gemmi/structure.hpp>"
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
