"""
coot_types.py — find Coot type definitions in the source tree.

Companion to call_sites.py: when a function takes/returns a `coot::Foo`, the
LLM has no idea what `Foo`'s fields are, how to construct one, or what
methods it offers. This module pulls the verbatim struct/class/enum/typedef
definition out of the header where it lives, namespace-aware, and lets
prompts.py inject it.

There is no LLM in this pipeline. Same trick as `call_sites`: ground the
model in real, in-tree code that it can mirror or read fields off.

Public API
----------
    find_coot_types(source_code)              -> list[CootType]
    format_coot_types_for_prompt(source_code) -> str
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

import report


# ── Tunables ─────────────────────────────────────────────────────────────────

_HEADER_EXTS = (".hh", ".h", ".hpp", ".hxx")

# Same exclusions as call_sites: skip generated dirs, this tooling tree, and
# external build outputs.
_EXCLUDE_DIR_NAMES = {
    ".git", "node_modules", ".venv", "venv", "__pycache__",
    "mmdb-recon", "mmdb-refactor-ui", "macos-build",
}
_EXCLUDE_DIR_PREFIXES = ("build",)

# Per-type body cap. Some coot structs are huge (`molecule_class_info_t` is
# thousands of lines). Truncate any single body that exceeds this so one fat
# type can't blow the prompt window.
_MAX_BODY_LINES = 80

# Total cap for the whole "Coot type definitions" section. We emit types in
# the order they appear in the function source (so signature types come
# first); once the budget is gone, we stop and emit a marker.
_MAX_TOTAL_LINES = 400

_DEFAULT_MAX_TYPES = 8


# ── Data type ────────────────────────────────────────────────────────────────

@dataclass
class CootType:
    fqn: str        # fully-qualified name, e.g. "coot::minimol::residue"
    kind: str       # "class" | "struct" | "enum" | "typedef"
    rel_path: str
    lineno: int
    body: str       # raw text from the header (full body for class/struct/enum,
                    # single declaration line for typedef/using)


# ── Module-level caches ──────────────────────────────────────────────────────

_file_list_cache: list[Path] | None = None
_index: dict[str, CootType] | None = None


# ── Comment / string scrubbing ───────────────────────────────────────────────
#
# We rewrite comments and string/char literals into spaces (preserving every
# original offset and newline) BEFORE running the structural regex. This stops
# `// class Foo` and `"namespace bar {"` style content from tricking the
# scanner into thinking there's a real declaration. Offsets are preserved so
# we can still slice the original `text` to recover bodies.

_LINE_COMMENT_RE  = re.compile(r"//[^\n]*")
_BLOCK_COMMENT_RE = re.compile(r"/\*[\s\S]*?\*/")
_STRING_LIT_RE    = re.compile(r'"(?:\\.|[^"\\\n])*"')
_CHAR_LIT_RE      = re.compile(r"'(?:\\.|[^'\\\n])*'")


def _blank(m: re.Match) -> str:
    s = m.group()
    return "".join(" " if c != "\n" else "\n" for c in s)


def _scrub(text: str) -> str:
    text = _BLOCK_COMMENT_RE.sub(_blank, text)
    text = _LINE_COMMENT_RE.sub(_blank, text)
    text = _STRING_LIT_RE.sub(_blank, text)
    text = _CHAR_LIT_RE.sub(_blank, text)
    return text


# ── Token regex ──────────────────────────────────────────────────────────────
#
# A single combined alternation finds every "interesting" token in one pass:
# namespaces, class/struct/enum bodies, typedef/using aliases, and bare
# braces (so we can track depth). Each alternative is a named group so we can
# dispatch on `m.group('group_name')` cheaply.

_TOKEN_RE = re.compile(
    # `namespace coot {` and C++17 `namespace coot::minimol {`
    r"(?P<ns_def>\bnamespace\s+(?P<ns_name>\w+(?:\s*::\s*\w+)*)\s*\{)"
    r"|(?P<ns_anon>\bnamespace\s*\{)"
    # `class Foo : public Bar {` / `struct Foo final {`
    r"|(?P<class_def>\b(?P<class_kind>class|struct)\s+(?P<class_name>\w+)\b[^;{}]*?\{)"
    # `enum Foo {` / `enum class Foo : int {`
    r"|(?P<enum_def>\benum\s+(?:class\s+|struct\s+)?(?P<enum_name>\w+)\b[^;{}]*?\{)"
    # `typedef T name;`
    r"|(?P<typedef_def>\btypedef\s+(?P<typedef_target>[^;{}]+?)\s+(?P<typedef_name>\w+)\s*(?:\[[^\]]*\])?\s*;)"
    # `using name = T;`
    r"|(?P<using_def>\busing\s+(?P<using_name>\w+)\s*=\s*(?P<using_target>[^;{}]+);)"
    # bare braces — only counted when none of the above matched
    r"|(?P<obrace>\{)"
    r"|(?P<cbrace>\})"
)


def _find_matching_close(cleaned: str, start: int) -> int | None:
    """Walk forward from `start` (just past an open brace) and return the
    offset of the matching close brace, or None if unbalanced."""
    depth = 1
    n = len(cleaned)
    i = start
    while i < n:
        ch = cleaned[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return None


def _qualify(scope_stack: list[tuple[str, int]], name: str) -> str:
    """Build a fully-qualified name from a scope stack + leaf name."""
    parts: list[str] = []
    for s_name, _ in scope_stack:
        if s_name:
            # `namespace coot::minimol {` lives as one entry but expands to
            # two FQN components.
            parts.extend(p.strip() for p in s_name.split("::"))
    parts.append(name)
    return "::".join(parts)


def _scan_text(text: str, rel_path: str) -> list[CootType]:
    """Find every namespace-qualified type defined in this header text."""
    cleaned = _scrub(text)
    out: list[CootType] = []

    # `(name_or_empty, depth_inside_scope)` for both namespaces and classes.
    # Classes go on the stack so nested types pick up `Outer::Inner`.
    scope_stack: list[tuple[str, int]] = []
    depth = 0

    for m in _TOKEN_RE.finditer(cleaned):
        if m.group("obrace"):
            depth += 1
            continue
        if m.group("cbrace"):
            depth -= 1
            while scope_stack and scope_stack[-1][1] > depth:
                scope_stack.pop()
            continue

        if m.group("ns_def"):
            depth += 1
            scope_stack.append((m.group("ns_name"), depth))
            continue
        if m.group("ns_anon"):
            depth += 1
            scope_stack.append(("", depth))
            continue

        if m.group("class_def"):
            depth += 1
            name = m.group("class_name")
            kind = m.group("class_kind")
            close = _find_matching_close(cleaned, m.end())
            if close is not None:
                body   = text[m.start():close + 2]   # include `};`
                lineno = text.count("\n", 0, m.start()) + 1
                out.append(CootType(
                    fqn=_qualify(scope_stack, name),
                    kind=kind,
                    rel_path=rel_path,
                    lineno=lineno,
                    body=body,
                ))
            # Push so any nested types we encounter while iterating get the
            # correct `Outer::Inner` qualified name.
            scope_stack.append((name, depth))
            continue

        if m.group("enum_def"):
            depth += 1
            name  = m.group("enum_name")
            close = _find_matching_close(cleaned, m.end())
            if close is not None:
                body   = text[m.start():close + 2]
                lineno = text.count("\n", 0, m.start()) + 1
                out.append(CootType(
                    fqn=_qualify(scope_stack, name),
                    kind="enum",
                    rel_path=rel_path,
                    lineno=lineno,
                    body=body,
                ))
            scope_stack.append((name, depth))
            continue

        if m.group("typedef_def"):
            name = m.group("typedef_name")
            out.append(CootType(
                fqn=_qualify(scope_stack, name),
                kind="typedef",
                rel_path=rel_path,
                lineno=text.count("\n", 0, m.start()) + 1,
                body=text[m.start():m.end()],
            ))
            continue

        if m.group("using_def"):
            name = m.group("using_name")
            out.append(CootType(
                fqn=_qualify(scope_stack, name),
                kind="typedef",
                rel_path=rel_path,
                lineno=text.count("\n", 0, m.start()) + 1,
                body=text[m.start():m.end()],
            ))
            continue

    return out


# ── Index build ──────────────────────────────────────────────────────────────

def _gather_header_files() -> list[Path]:
    global _file_list_cache
    if _file_list_cache is not None:
        return _file_list_cache
    root = report._coot_root
    if not root:
        _file_list_cache = []
        return _file_list_cache
    files: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [
            d for d in dirnames
            if d not in _EXCLUDE_DIR_NAMES
            and not d.startswith(".")
            and not any(d.startswith(p) for p in _EXCLUDE_DIR_PREFIXES)
        ]
        for fn in filenames:
            if fn.endswith(_HEADER_EXTS):
                files.append(Path(dirpath) / fn)
    _file_list_cache = files
    return files


def _build_index() -> dict[str, CootType]:
    files = _gather_header_files()
    if not files:
        return {}
    root = Path(report._coot_root)
    index: dict[str, CootType] = {}
    for path in files:
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        try:
            rel = str(path.relative_to(root))
        except ValueError:
            rel = str(path)
        for t in _scan_text(text, rel):
            # Only keep `coot::*` types — everything else is irrelevant to
            # functions in the coot tree (and we'd otherwise capture every
            # struct in every system header we walked into).
            if t.fqn.startswith("coot::"):
                index.setdefault(t.fqn, t)
    print(f"coot_types: indexed {len(index)} coot types from {len(files)} headers")
    return index


def _get_index() -> dict[str, CootType]:
    global _index
    if _index is None:
        _index = _build_index()
    return _index


# ── Public API ───────────────────────────────────────────────────────────────

# Match `coot::Foo`, `coot::Foo::Bar`, `coot::ns::Type::method` etc. The
# top-level pattern captures everything after `coot::`; we then break it
# down into prefixes when looking up.
_COOT_REF_RE = re.compile(r"\bcoot::([A-Za-z_][A-Za-z0-9_]*(?:::[A-Za-z_][A-Za-z0-9_]*)*)\b")


def find_coot_types(
    source_code: str,
    max_types: int = _DEFAULT_MAX_TYPES,
) -> list[CootType]:
    """Return verbatim coot type definitions referenced by `source_code`.

    The function source is scanned for `coot::*` identifiers in document order
    (so signature types come first). For each identifier we try the longest
    `coot::a::b::c` prefix that resolves to a type in the index, then shorter
    prefixes — that handles `coot::residue_spec_t::null_residue` (where the
    type is the prefix and `null_residue` is a member).
    """
    if not source_code:
        return []
    index = _get_index()
    if not index:
        return []

    seen: set[str] = set()
    out: list[CootType] = []
    for m in _COOT_REF_RE.finditer(source_code):
        path  = m.group(1)
        parts = path.split("::")
        for n in range(len(parts), 0, -1):
            fqn = "coot::" + "::".join(parts[:n])
            if fqn in seen:
                break
            if fqn in index:
                seen.add(fqn)
                out.append(index[fqn])
                if len(out) >= max_types:
                    return out
                break
    return out


def _truncate_body(body: str) -> str:
    lines = body.splitlines()
    if len(lines) <= _MAX_BODY_LINES:
        return body
    head  = lines[: _MAX_BODY_LINES - 1]
    extra = len(lines) - len(head)
    return "\n".join(head) + f"\n// ... ({extra} more lines elided)"


def format_coot_types_for_prompt(
    source_code: str,
    max_types: int = _DEFAULT_MAX_TYPES,
) -> str:
    """Render the matched coot types as a markdown prompt section, capped to
    `_MAX_TOTAL_LINES` lines total. Returns '' if nothing resolved."""
    types = find_coot_types(source_code, max_types=max_types)
    if not types:
        return ""

    blocks: list[str] = []
    total_lines = 0
    for i, t in enumerate(types):
        body = _truncate_body(t.body)
        block_lines = body.count("\n") + 1
        if total_lines + block_lines > _MAX_TOTAL_LINES and blocks:
            blocks.append(
                f"// ... ({len(types) - i} more types omitted to fit prompt budget)"
            )
            break
        total_lines += block_lines
        blocks.append(
            f"### `{t.fqn}` ({t.kind})\n"
            f"// {t.rel_path}:{t.lineno}\n"
            f"```cpp\n{body}\n```"
        )

    body_md = "\n\n".join(blocks)
    return (
        "\n\n## Coot type definitions used by this function\n\n"
        "The function references the following Coot internal types. Their "
        "verbatim definitions are pulled from the coot source tree so you "
        "know exactly what fields they have, what constructors they accept, "
        "and which methods you can call. Use ONLY fields and methods that "
        "appear in these definitions — do not invent members. If a field is "
        "named `chain_id`, it is `chain_id`, not `chainID` or `getChainID()`.\n\n"
        f"{body_md}\n"
    )
