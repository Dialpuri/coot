"""
doc_extractor.py — Extract relevant class/method sections from LLM-generated
markdown docs (produced by document_methods.py) given a list of MMDB symbols.

Markdown structure assumed:
  ### `file_stem`          ← file group (level 3)
  ## `ClassName`           ← class / struct (level 2)
  ### Methods
  #### `return method()`   ← individual method (level 4)
  ## Free Functions
  #### `return fn()`

Real MMDB symbol formats seen in the AST report:
  "class mmdb::Chain"                      → class Chain
  "mmdb::Chain"                            → type Chain
  "mmdb::Chain *"                          → pointer to Chain
  "mmdb::Chain::GetNumberOfResidues"       → method on Chain
  "mmdb::PPAtom"                           → pointer-pointer type → Atom
  "mmdb::mat44"                            → plain type

Public API
----------
  load_docs(path)                        -> str
  extract_for_symbols(markdown, symbols) -> str
"""

from __future__ import annotations

import re
from pathlib import Path


# ---------------------------------------------------------------------------
# Markdown parsing
# ---------------------------------------------------------------------------

_H2 = re.compile(r"^## (.+)$", re.MULTILINE)
_H4 = re.compile(r"^#### (.+)$", re.MULTILINE)


def _split_by_h2(markdown: str) -> list[tuple[str, str]]:
    """[(heading_text, section_body), ...] for every ## section."""
    sections: list[tuple[str, str]] = []
    matches = list(_H2.finditer(markdown))
    for i, m in enumerate(matches):
        heading = m.group(1).strip().strip("`")
        start   = m.end()
        end     = matches[i + 1].start() if i + 1 < len(matches) else len(markdown)
        sections.append((heading, markdown[start:end]))
    return sections


def _split_h4_blocks(body: str) -> list[tuple[str, str]]:
    """[(method_heading, full_block_text), ...] for every #### inside a class body."""
    blocks: list[tuple[str, str]] = []
    matches = list(_H4.finditer(body))
    for i, m in enumerate(matches):
        heading = m.group(1).strip().strip("`")
        start   = m.start()
        end     = matches[i + 1].start() if i + 1 < len(matches) else len(body)
        blocks.append((heading, body[start:end]))
    return blocks


def _method_name_from_heading(heading: str) -> str:
    """Extract the bare method name from '#### int Foo::Bar(int x)'."""
    m = re.search(r"(\w+)\s*\(", heading)
    return m.group(1) if m else heading


# ---------------------------------------------------------------------------
# Symbol normalisation
# ---------------------------------------------------------------------------

def _parse_symbol(raw: str) -> tuple[str | None, str | None]:
    """
    Parse a raw MMDB symbol string into (class_name, method_name).
    Either component may be None.

    Examples
    --------
    "class mmdb::Chain"              -> ("Chain", None)
    "mmdb::Chain"                    -> ("Chain", None)
    "mmdb::Chain *"                  -> ("Chain", None)
    "mmdb::Chain::GetNumberOfResidues" -> ("Chain", "GetNumberOfResidues")
    "mmdb::PPAtom"                   -> ("Atom",  None)   [PP-prefix stripped]
    "mmdb::PAtom"                    -> ("Atom",  None)   [P-prefix stripped]
    "mmdb::mat44"                    -> ("mat44", None)
    """
    s = raw.strip()

    # Strip 'class ' / 'struct ' prefix
    s = re.sub(r"^(class|struct)\s+", "", s)

    # Strip trailing pointer markers: ' *', '*', ' &', '&'
    s = re.sub(r"\s*[*&]+\s*$", "", s).strip()

    # Split on '::'
    parts = [p for p in s.split("::") if p]  # drop empty strings from leading ::

    # Discard 'mmdb' / 'mmdb2' namespace token
    parts = [p for p in parts if p.lower() not in ("mmdb", "mmdb2")]

    if not parts:
        return None, None

    if len(parts) == 1:
        # Single token → could be a type (Atom, mat44) or pointer type (PPAtom)
        bare = parts[0]
        # Strip P / PP pointer-type prefixes common in MMDB
        stripped = re.sub(r"^P{1,2}(?=[A-Z])", "", bare)
        class_name = stripped if stripped != bare else bare
        return class_name, None

    # Two or more tokens → last token is likely a method, rest is class path
    method_name = parts[-1]
    class_name  = parts[-2]  # e.g. Chain from mmdb::Chain::GetNumberOfResidues

    # Strip P/PP prefix from class name too just in case
    class_name = re.sub(r"^P{1,2}(?=[A-Z])", "", class_name)

    return class_name, method_name


def _heading_matches_class(heading: str, class_name: str) -> bool:
    """
    True if the section heading refers to the given class.

    Heading examples: 'mmdb::Atom', 'mmdb::Chain', 'Matrix', 'Residue'
    We match the bare class name (case-insensitive) against the last component
    of the heading.
    """
    # Strip backticks
    h = heading.strip("`").strip()
    # Get the last :: component of the heading (e.g. 'mmdb::Atom' -> 'Atom')
    h_bare = h.split("::")[-1].strip()
    return h_bare.lower() == class_name.lower()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_docs(path: str | Path) -> str:
    """Read the markdown file and return its contents, or '' if missing."""
    p = Path(path)
    if not p.exists():
        return ""
    return p.read_text(errors="replace")


def extract_for_symbols(markdown: str, symbols: list[str]) -> str:
    """
    Given the full docs markdown and a list of raw MMDB symbol strings,
    return a trimmed markdown string with only the relevant class sections
    (and, where a specific method is named, only that method's block).

    Algorithm
    ---------
    1. Parse each symbol into (class_name, method_name).
    2. Group by class_name: build a map  class_name -> set of wanted methods
       (empty set = want the whole class section).
    3. For each ## heading in the markdown, check if it matches a wanted class.
       - If matched and no specific methods wanted → include full section.
       - If matched and specific methods wanted → include preamble + matching
         method blocks only. Fall back to full section if no methods matched.
    """
    if not markdown or not symbols:
        return ""

    # Build wanted map: class_name (lower) -> set of method names (lower)
    # An empty set means "want the whole class section"
    wanted: dict[str, set[str]] = {}
    for raw in symbols:
        class_name, method_name = _parse_symbol(raw)
        if not class_name:
            continue
        key = class_name.lower()
        if key not in wanted:
            wanted[key] = set()
        if method_name:
            wanted[key].add(method_name.lower())
        # If we later encounter the same class without a method (i.e. want whole
        # class), mark with None sentinel by clearing the set
        # (handled below — an empty set after all parsing means "whole class")

    # Second pass: if a class appears both with and without a method,
    # want the full class (empty set = full class).
    for raw in symbols:
        class_name, method_name = _parse_symbol(raw)
        if class_name and method_name is None:
            wanted[class_name.lower()] = set()  # override: want full section

    if not wanted:
        return ""

    sections = _split_by_h2(markdown)
    collected: list[str] = []
    seen: set[str] = set()

    for heading, body in sections:
        h_bare = heading.strip("`").split("::")[-1].strip().lower()
        if h_bare not in wanted:
            continue
        if heading in seen:
            continue
        seen.add(heading)

        wanted_methods = wanted[h_bare]

        if not wanted_methods:
            # Include the entire class section
            collected.append(f"## `{heading}`\n{body}".rstrip())
            continue

        # Include only the matching method blocks
        method_blocks = _split_h4_blocks(body)
        matched = [
            block
            for mh, block in method_blocks
            if _method_name_from_heading(mh).lower() in wanted_methods
        ]

        if matched:
            first_h4 = _H4.search(body)
            preamble  = body[:first_h4.start()] if first_h4 else ""
            collected.append(
                (f"## `{heading}`\n{preamble}" + "".join(matched)).rstrip()
            )
        else:
            # Method not found in this section — include the full section
            # so the LLM at least has the class context
            collected.append(f"## `{heading}`\n{body}".rstrip())

    return "\n\n---\n\n".join(collected)
