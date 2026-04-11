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
from call_sites import format_call_sites_for_prompt
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

# System prompt for MMDB-only test generation. The default STRATEGY_CONTEXT is
# all about migrating MMDB → Gemmi, so it pushes the model to set up a gemmi
# structure even when the user explicitly asked for an MMDB test. This prompt
# pins the model to MMDB-only output.
MMDB_TEST_SYSTEM_CONTEXT: str = (
    "You are an expert C++ developer writing a Google Test for an EXISTING "
    "MMDB2-based function. This is NOT a refactor and NOT a migration.\n\n"
    "CRITICAL RULES:\n"
    "- Use MMDB2 types ONLY. Do NOT emit `gemmi::` types, "
    "`#include <gemmi/...>`, `gemmi::read_structure_file`, or any other "
    "Gemmi code anywhere in the test. The test must be pure MMDB.\n"
    "- The ONLY MMDB header you need is `#include <mmdb2/mmdb_manager.h>` "
    "plus the function's own header. Do not include `mmdb_model.h`, "
    "`mmdb_chain.h`, `mmdb_atom.h`, `mmdb_root.h`, or any other mmdb "
    "header — those paths do not exist as standalone public headers.\n"
    "- Every MMDB type MUST be qualified with `mmdb::`: `mmdb::Manager`, "
    "`mmdb::Model*`, `mmdb::Chain*`, `mmdb::Residue*`, `mmdb::Atom*`, "
    "`mmdb::Error_NoError`. There is no `using namespace mmdb;` — "
    "unqualified types will not compile.\n"
    "- Load molecular data by having `mmdb::Manager::ReadCoorFile` parse a "
    "real PDB file. Do NOT hand-build synthetic mmdb structures.\n"
    "- Do NOT mock the function under test. Call its real declaration "
    "directly.\n"
)

# System prompt for the "both" target — produces both an MMDB test and a
# Gemmi test in one shot. Here gemmi IS allowed (in the gemmi block only),
# so we keep the migration-style framing.
BOTH_TEST_SYSTEM_CONTEXT: str = (
    "You are an expert C++ developer writing a pair of Google Tests for a "
    "function that exists in MMDB2 today and will be ported to Gemmi.\n\n"
    "Output two TEST() blocks:\n"
    "  1. An MMDB block that uses MMDB2 types ONLY (no gemmi anywhere).\n"
    "  2. A Gemmi block that uses Gemmi types ONLY (no mmdb anywhere).\n\n"
    "Do not mix the two libraries inside a single TEST() block. Both blocks "
    "load a real PDB file from disk and assert against the measured oracle "
    "values supplied in the prompt — never invent expected numbers.\n"
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


def _call_sites_context(symbols: list[str]) -> str:
    """Wrap `call_sites.format_call_sites_for_prompt` so a missing index or
    walk failure can never break prompt construction."""
    try:
        return format_call_sites_for_prompt(symbols)
    except Exception as ex:
        print(f"call_sites: failed to gather examples ({ex}) — continuing without them")
        return ""


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


# ── Source-file include scraping ──────────────────────────────────────────────
#
# Rather than relying on the AST tool to emit captured #include directives,
# we just read the .cc file off disk and parse its include lines directly.
# This is simpler, always up-to-date with the source tree, and needs no
# rebuild/rescan cycle.

_INCLUDE_RE = re.compile(r'^\s*#\s*include\s*([<"][^>"]+[>"])')


def read_source_includes(rel_source_path: str) -> list[str]:
    """Return the list of #include directives (in order, deduplicated) that
    appear at the top of the given source file, as raw spellings like
    `<mmdb2/mmdb_manager.h>` or `"coot-coord-utils.hh"`.

    Stops scanning as soon as we see a non-include, non-blank, non-comment
    line that looks like real code — we only want the header block.
    """
    import report
    if not rel_source_path:
        return []
    src = Path(report._coot_root) / rel_source_path
    if not src.exists():
        return []
    seen: set[str] = set()
    out: list[str] = []
    try:
        with src.open("r", errors="replace") as fh:
            in_block_comment = False
            for raw in fh:
                line = raw.rstrip("\n")
                stripped = line.strip()
                # Handle /* ... */ block comments that span lines
                if in_block_comment:
                    if "*/" in stripped:
                        in_block_comment = False
                    continue
                if stripped.startswith("/*") and "*/" not in stripped:
                    in_block_comment = True
                    continue
                if not stripped:
                    continue
                if stripped.startswith("//") or stripped.startswith("/*"):
                    continue
                # Preprocessor directives are fine to skip past
                m = _INCLUDE_RE.match(line)
                if m:
                    spelling = m.group(1)
                    if spelling not in seen:
                        seen.add(spelling)
                        out.append(spelling)
                    continue
                if stripped.startswith("#"):
                    # other preprocessor (define, if, pragma, etc.) — keep going
                    continue
                # First line of real code → include block is over
                break
    except OSError:
        return []
    return out


def _format_file_includes(rel_source_path: str) -> str:
    """Render the source file's real #include block as a prompt section.

    This tells the LLM the exact headers the original translation unit uses,
    so it doesn't invent include paths that don't exist on disk.
    """
    return ""

    includes = read_source_includes(rel_source_path)
    if not includes:
        return ""
    block = "\n".join(f"#include {inc}" for inc in includes)
    return (
        "\n\n## Real includes from the source file\n\n"
        f"The translation unit `{rel_source_path}` begins with these "
        "`#include` directives. Use EXACTLY these spellings in your probe — "
        "do not invent header paths, do not drop the leading directory, do "
        "not switch `<...>` to `\"...\"` or vice versa. Copy the ones you "
        "actually need verbatim:\n\n"
        f"```cpp\n{block}\n```\n"
    )


def _format_header_include(rel_source_path: str) -> str:
    """Render a '## Where the function lives' prompt section, or empty string.

    Tells the LLM explicitly which header to `#include` and stresses that the
    function is real — not to be mocked. Also pins down the MMDB header so
    the LLM stops dragging in `mmdb_model.h`, `mmdb_chain.h`, etc., none of
    which exist as standalone public headers.
    """
    mmdb_rule = (
        "The ONLY MMDB header you need is `#include <mmdb2/mmdb_manager.h>`. "
        "It declares `mmdb::Manager`, `mmdb::Model`, `mmdb::Chain`, "
        "`mmdb::Residue`, `mmdb::Atom`, and the file-loading API. Do NOT "
        "`#include <mmdb2/mmdb_model.h>`, `<mmdb2/mmdb_chain.h>`, "
        "`<mmdb2/mmdb_atom.h>`, `<mmdb2/mmdb_root.h>`, or any other mmdb "
        "header — those paths do not exist.\n\n"
        "Every MMDB type MUST be written with the `mmdb::` namespace prefix "
        "on first use and every use after that — write `mmdb::Manager`, "
        "`mmdb::Model*`, `mmdb::Chain*`, `mmdb::Residue*`, `mmdb::Atom*`, "
        "`mmdb::Error_NoError`. There is no `using namespace mmdb;` in coot "
        "test code, so an unqualified `Manager` or `Chain*` will not "
        "compile. Constants like `mmdb::Error_NoError` and `mmdb::PPAtom` "
        "are also in the `mmdb` namespace."
    )

    if not rel_source_path:
        return (
            "\n\n## Where the function lives\n\n"
            f"The function is defined in the coot tree. Do NOT mock it, "
            "reimplement it, or declare a local stub — `#include` its real "
            "header and call it directly.\n\n"
            f"{mmdb_rule}\n"
        )
    header = find_header_for_source(rel_source_path)
    if not header:
        return (
            "\n\n## Where the function lives\n\n"
            f"The function is defined in `{rel_source_path}`. Do NOT mock it, "
            "reimplement it, or declare a local stub — `#include` the real "
            "declaration from its header (in the same directory as the .cc "
            "file) and call it directly. The `-I<coot root>` include path is "
            "already on the compile command.\n\n"
            f"{mmdb_rule}\n"
        )
    return (
        "\n\n## Where the function lives\n\n"
        f"The function is declared in `{header}` (defined in `{rel_source_path}`).\n\n"
        f"You MUST include this header and call the real function — do NOT mock it, "
        "reimplement it, or declare a local stub:\n\n"
        f"```cpp\n#include \"{header}\"\n```\n\n"
        "The coot source tree is already on the compile command's include path "
        f"(`-I<coot root>`), so this include will resolve without any extra flags.\n\n"
        f"{mmdb_rule}\n"
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


def system_context_for_test_target(target: str) -> str:
    """Return the right system prompt for an `/api/generate-test` target.

    - "mmdb"  → pure MMDB framing, no gemmi at all
    - "gemmi" → MMDB→Gemmi strategy context (the gemmi test is the migration)
    - "both"  → dual-block framing where each block stays in its own library
    """
    if target == "mmdb":
        return MMDB_TEST_SYSTEM_CONTEXT
    if target == "both":
        return BOTH_TEST_SYSTEM_CONTEXT
    return STRATEGY_CONTEXT


def build_generate_test(function_name: str, source_code: str, mmdb_symbols: list[str],
                        target: str, additional_instructions: str = "",
                        oracle_output: str = "",
                        rel_source_path: str = "",
                        probe_source: str = "",
                        probe_pdb_path: str = "") -> str:
    """Build a test-generation prompt for target='mmdb'|'gemmi'|'both'.

    Inputs from a successful oracle run
    -----------------------------------
    `oracle_output`  — captured PROBE: lines from the real run, used as the
                       expected values inside `EXPECT_*` assertions.
    `probe_source`   — the C++ source of the probe `main()` that produced
                       those values. The MMDB test should be that source
                       wrapped in a `TEST()` block, not a freshly-imagined
                       loader. Including it stops the LLM from inventing
                       its own (broken) data-loading code.
    `probe_pdb_path` — the PDB file the probe used; the test must read the
                       same file so its assertions match.

    If `rel_source_path` is provided, we look up the matching header and tell
    the LLM to `#include` it so it calls the real function instead of mocking.
    """
    common = dict(
        function_name=function_name,
        symbols_list=", ".join(mmdb_symbols) if mmdb_symbols else "none",
        additional_instructions=_additional(additional_instructions),
        api_context=_api_context(mmdb_symbols),
        call_sites=_call_sites_context(mmdb_symbols),
        source_code=source_code,
        gtest_style=GTEST_STYLE,
        oracle_values=_format_oracle(oracle_output),
        header_include=_format_header_include(rel_source_path),
        file_includes=_format_file_includes(rel_source_path),
        probe_reference=_format_probe_reference(probe_source, probe_pdb_path),
    )
    if target == "mmdb":
        return _T_TEST_MMDB.substitute(**common)
    if target == "gemmi":
        return _T_TEST_GEMMI.substitute(**common)
    return _T_TEST_MMDB.substitute(**common)


# Symbols the probe ALWAYS needs docs for, regardless of what the function
# under test happens to touch. The probe has to open a PDB file and walk the
# mmdb hierarchy down to whatever object the function expects — but the
# function itself usually receives its arguments already-constructed, so its
# own AST-derived `mmdb_symbols` rarely mention `ReadCoorFile`, `GetModel`,
# `GetChain`, etc. We force those docs into the probe's API context so the
# LLM actually knows how to load a file and navigate the hierarchy.
#
# These names are matched against the doc headings by doc_extractor; anything
# that doesn't resolve is silently dropped, so over-listing is safe.
_PROBE_ESSENTIAL_SYMBOLS: list[str] = [
    # File I/O (lives on Root / Manager)
    "mmdb::Manager",
    "mmdb::Root::ReadCoorFile",
    "mmdb::Root::ReadPDBASCII",
    "mmdb::Root::ReadCIFASCII",
    "mmdb::Root::GetErrorDescription",
    # Hierarchy navigation: Manager/CoorManager → Model → Chain → Residue → Atom
    "mmdb::CoorManager::GetNumberOfModels",
    "mmdb::CoorManager::GetModel",
    "mmdb::CoorManager::GetFirstDefaultModel",
    "mmdb::Model::GetNumberOfChains",
    "mmdb::Model::GetChain",
    "mmdb::Chain::GetChainID",
    "mmdb::Chain::GetNumberOfResidues",
    "mmdb::Chain::GetResidue",
    "mmdb::Residue::GetNumberOfAtoms",
    "mmdb::Residue::GetAtom",
    "mmdb::Residue::GetResName",
    "mmdb::Residue::GetSeqNum",
    "mmdb::Atom::GetAtomName",
    "mmdb::Atom::GetChainID",
    "mmdb::Atom::GetResidue",
]


def build_probe_mmdb(function_name: str, source_code: str,
                     mmdb_symbols: list[str], additional_instructions: str = "",
                     rel_source_path: str = "", pdb_path: str = "") -> str:
    """Build a prompt asking the LLM to write an MMDB probe main().

    The probe needs to load a PDB file and walk the mmdb hierarchy down to
    the right object, so we merge `_PROBE_ESSENTIAL_SYMBOLS` into the
    caller-supplied symbol list before resolving docs. This ensures the LLM
    always sees the `ReadCoorFile` / `GetModel` / `GetChain` / etc. API
    reference, even if none of those appear in the function's own AST symbols.
    """
    merged_symbols: list[str] = list(mmdb_symbols)
    seen = set(merged_symbols)
    for sym in _PROBE_ESSENTIAL_SYMBOLS:
        if sym not in seen:
            merged_symbols.append(sym)
            seen.add(sym)

    return _T_PROBE_MMDB.substitute(
        function_name=function_name,
        symbols_list=", ".join(mmdb_symbols) if mmdb_symbols else "none",
        additional_instructions=_additional(additional_instructions),
        api_context=_api_context(merged_symbols),
        call_sites=_call_sites_context(merged_symbols),
        source_code=source_code,
        header_include=_format_header_include(rel_source_path),
        file_includes=_format_file_includes(rel_source_path),
        pdb_path=pdb_path,
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


def _format_probe_reference(probe_source: str, pdb_path: str) -> str:
    """Render the working probe `main()` as a reference for the test prompt.

    The probe has already been compiled and run successfully — it loaded the
    PDB, walked the mmdb hierarchy to the right object, and called the
    function under test. The MMDB test should reuse this exact loading code
    inside a TEST() block instead of inventing its own. Showing the LLM the
    working source has been the most reliable way to stop it from mocking
    things or guessing the API.
    """
    src = (probe_source or "").strip()
    if not src:
        return ""
    pdb_note = (
        f"The probe loaded `{pdb_path}`. Your MMDB test must read the SAME "
        "file so its assertions match the measured values shown above.\n\n"
        if pdb_path else ""
    )
    return (
        "\n\n## Working probe (already compiled and run successfully)\n\n"
        "The C++ source below is a `main()` that loaded a real PDB file, "
        "walked the mmdb hierarchy, and called the function under test to "
        "produce the measured values listed above. Treat it as ground truth.\n\n"
        f"{pdb_note}"
        "When writing the MMDB test:\n"
        "  - COPY the data-loading and function-call lines from this probe "
        "into your `TEST()` body verbatim. Do not invent a different loader.\n"
        "  - Drop the `PROBE:` `std::cout` lines and replace them with "
        "`EXPECT_EQ` / `EXPECT_NEAR` / `EXPECT_STREQ` against the measured "
        "values.\n"
        "  - Use ONLY the headers the probe uses — typically just "
        "`#include <mmdb2/mmdb_manager.h>` plus the function's own header.\n\n"
        f"```cpp\n{src}\n```\n"
    )


def build_fix(fn_name: str, variant: str, current_code: str,
              error: str, attempt: int) -> str:
    if variant == "mmdb":
        inc_hint = (
            "Required headers (use ONLY these — do not add other mmdb headers):\n"
            "#include <gtest/gtest.h>\n"
            "#include <mmdb2/mmdb_manager.h>\n"
            "\n"
            "Every MMDB type MUST be qualified with `mmdb::` — write "
            "`mmdb::Manager`, `mmdb::Model*`, `mmdb::Chain*`, "
            "`mmdb::Residue*`, `mmdb::Atom*`, `mmdb::Error_NoError`. "
            "Unqualified `Manager` / `Chain*` will not compile."
        )
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
