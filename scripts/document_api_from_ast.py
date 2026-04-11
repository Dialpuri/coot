#!/usr/bin/env python3
"""
Generate targeted MMDB API documentation driven by the AST usage report.

Instead of documenting every method in every file, this script works backwards:
  1. Reads the AST report to find every Class::Method pair actually used in the
     codebase being refactored.
  2. Maps each class to its MMDB2 header file and extracts the class declaration.
  3. Sends the header (or a filtered subset) + a list of the used methods to the
     LLM, asking it to document *only* those methods.
  4. Assembles a single markdown file compatible with doc_extractor.py.

This catches methods that document_methods.py misses because they are:
  - Only declared in headers (no paired .cpp)
  - Inherited from abstract base classes (ProModel, Root, SelManager, etc.)
  - Defined in header-only or template files

Usage:
    # Anthropic Claude (default)
    python document_api_from_ast.py \\
        --report mmdb-recon/mmdb_usage_report.json \\
        --mmdb-include /opt/homebrew/Cellar/mmdb2/2.0.22/include/mmdb2 \\
        --output mmdb-recon/llm/mmdb_methods.md

    # Ollama
    python document_api_from_ast.py \\
        --report mmdb-recon/mmdb_usage_report.json \\
        --mmdb-include /opt/homebrew/Cellar/mmdb2/2.0.22/include/mmdb2 \\
        --provider ollama --model codellama:34b
"""

import argparse
import json
import os
import re
import sys
import time
from abc import ABC, abstractmethod
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

DEFAULT_ANTHROPIC_MODEL = "claude-opus-4-6"
DEFAULT_OLLAMA_MODEL    = "gemma4"
DEFAULT_OLLAMA_URL      = "http://localhost:11434"
RATE_LIMIT_DELAY        = 1.0


# ---------------------------------------------------------------------------
# LLM backends (identical to summarize_codebase.py / document_methods.py)
# ---------------------------------------------------------------------------

class LLMBackend(ABC):
    @abstractmethod
    def complete(self, system: str, user: str, max_tokens: int) -> str: ...

    @property
    @abstractmethod
    def label(self) -> str: ...


class AnthropicBackend(LLMBackend):
    def __init__(self, model: str, api_key: str, think: bool = True) -> None:
        import anthropic
        self._client    = anthropic.Anthropic(api_key=api_key)
        self._model     = model
        self._anthropic = anthropic
        self._think     = think

    @property
    def label(self) -> str:
        suffix = "" if self._think else " (thinking off)"
        return f"Claude ({self._model}){suffix}"

    def complete(self, system: str, user: str, max_tokens: int) -> str:
        kwargs: dict = dict(
            model=self._model,
            max_tokens=max_tokens,
            system=system,
            messages=[{"role": "user", "content": user}],
        )
        if self._think:
            kwargs["thinking"] = {"type": "adaptive"}
        while True:
            try:
                with self._client.messages.stream(**kwargs) as stream:
                    final = stream.get_final_message()
                for block in reversed(final.content):
                    if block.type == "text":
                        return block.text
                return ""
            except self._anthropic.RateLimitError:
                print("rate limited — waiting 30s ...", end=" ", flush=True)
                time.sleep(30)


class OllamaBackend(LLMBackend):
    def __init__(self, model: str, base_url: str = DEFAULT_OLLAMA_URL,
                 think: bool = True) -> None:
        try:
            import httpx
        except ImportError:
            print("Error: httpx is required.  pip install httpx", file=sys.stderr)
            sys.exit(1)
        self._httpx    = httpx
        self._model    = model
        self._base_url = base_url.rstrip("/")
        self._think    = think
        self._check_connection()

    def _check_connection(self) -> None:
        try:
            self._httpx.get(f"{self._base_url}/api/tags", timeout=5).raise_for_status()
        except Exception as e:
            print(f"Error: Cannot reach Ollama at {self._base_url}\n  {e}", file=sys.stderr)
            sys.exit(1)

    @property
    def label(self) -> str:
        suffix = "" if self._think else " (thinking off)"
        return f"Ollama ({self._model} @ {self._base_url}){suffix}"

    def complete(self, system: str, user: str, max_tokens: int) -> str:
        # Doc generation is reasoning-heavy: classify each method as static vs
        # instance, infer semantics from name+signature. We enable `think: true`
        # by default so reasoning-capable models actually reason. Ollama silently
        # ignores `think` on models that don't support it, so this is safe.
        payload: dict = {
            "model": self._model,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user",   "content": user},
            ],
            "stream": True,
            "options": {"num_predict": max_tokens},
        }
        if self._think:
            payload["think"] = True

        chunks: list[str] = []
        with self._httpx.stream("POST", f"{self._base_url}/api/chat",
                                json=payload, timeout=600) as response:
            response.raise_for_status()
            for line in response.iter_lines():
                if not line:
                    continue
                try:
                    data = json.loads(line)
                except json.JSONDecodeError:
                    continue
                # Ollama puts the visible answer in message.content and the
                # reasoning trace (when `think: true`) in message.thinking.
                # We deliberately drop thinking — the doc generator only needs
                # the final markdown — but it still helps the model produce
                # better classifications.
                chunks.append(data.get("message", {}).get("content", ""))
                if data.get("done"):
                    break
        return "".join(chunks)


# ---------------------------------------------------------------------------
# Phase 1: Extract used symbols from AST report
# ---------------------------------------------------------------------------

def _normalise_symbol(raw: str) -> "tuple[str | None, str | None]":
    """
    Parse a raw AST symbol into (class_name, method_name).
    Returns (None, None) for unrecognised / irrelevant symbols.

    Input examples:
      "class mmdb::Chain"                   -> ("Chain", None)
      "mmdb::Chain *"                       -> ("Chain", None)
      "mmdb::Chain::GetNumberOfResidues"    -> ("Chain", "GetNumberOfResidues")
      "mmdb::PPAtom"                        -> ("Atom",  None)
      "vector<mmdb::Residue *>"             -> ignored (None, None)
      "map<int, mmdb::Atom *>"              -> ignored
    """
    s = raw.strip()

    # Ignore C++ STL containers and other compound types
    if re.search(r"[<>]", s):
        return None, None

    # Strip 'class '/'struct ' prefix
    s = re.sub(r"^(class|struct)\s+", "", s)
    # Strip trailing pointer / reference markers
    s = re.sub(r"\s*[*&]+\s*$", "", s).strip()

    # Split on '::'
    parts = [p for p in s.split("::") if p.lower() not in ("mmdb", "mmdb2", "")]
    if not parts:
        return None, None

    if len(parts) == 1:
        bare = parts[0]
        # Strip P/PP MMDB pointer-type prefixes
        stripped = re.sub(r"^P{1,2}(?=[A-Z])", "", bare)
        return (stripped if stripped != bare else bare), None

    # Two components: Class::Method
    class_name  = re.sub(r"^P{1,2}(?=[A-Z])", "", parts[-2])
    method_name = parts[-1]
    return class_name, method_name


def collect_used_symbols(report_path: str) -> dict[str, set[str]]:
    """
    Returns {class_name: {method_name, ...}} for every MMDB class/method
    referenced in the AST report.  An empty method set means the class is used
    as a type but no specific methods were identified.
    """
    report = json.loads(Path(report_path).read_text())
    used: dict[str, set[str]] = {}

    for file_entry in report["files"]:
        for fn in file_entry.get("functions", []):
            for raw_sym in fn.get("mmdb_symbols", []):
                cls, meth = _normalise_symbol(raw_sym)
                if not cls:
                    continue
                if cls not in used:
                    used[cls] = set()
                if meth:
                    used[cls].add(meth)

    return used


# ---------------------------------------------------------------------------
# Phase 2: Map class names to header files
# ---------------------------------------------------------------------------

def build_class_to_header_map(include_dir: str) -> dict[str, Path]:
    """Scan MMDB2 headers and return {class_name: header_path}."""
    mapping: dict[str, Path] = {}
    for hfile in Path(include_dir).glob("*.h"):
        content = hfile.read_text(errors="replace")
        for m in re.finditer(
            r"class\s+(?:MMDB_DL_EXPORT\s+)?(\w+)\s*[:{(]", content
        ):
            cls = m.group(1)
            # Keep first occurrence (most likely the primary definition)
            mapping.setdefault(cls, hfile)
    return mapping


# ---------------------------------------------------------------------------
# Phase 3: Extract class body from header
# ---------------------------------------------------------------------------

def extract_class_body(header_text: str, class_name: str) -> "str | None":
    """
    Extract the full class declaration from the header text using brace counting.
    Returns the raw C++ text from 'class ClassName' to the closing '};'.
    """
    pattern = rf"class\s+(?:MMDB_DL_EXPORT\s+)?{re.escape(class_name)}\s*[:({{]"
    m = re.search(pattern, header_text)
    if not m:
        return None

    # Find the opening brace from the match position
    try:
        brace_pos = header_text.index("{", m.start())
    except ValueError:
        return None

    depth = 0
    i = brace_pos
    while i < len(header_text):
        ch = header_text[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return header_text[m.start() : i + 2]  # include ';'
        i += 1
    return header_text[m.start():]  # unclosed — return remainder


def find_static_methods(class_body: str, wanted_methods: set[str]) -> set[str]:
    """Return the subset of `wanted_methods` declared `static` in this class.

    The model keeps writing `mmdb::Manager::GetModel(1)` for instance methods,
    so we need to know unambiguously which methods can legitimately be called
    via `ClassName::method(...)` and which must go through a receiver.

    Heuristic: for each method, look for a token sequence
        static  <stuff that is not ; { }>  methodName (
    inside the class body. The forbidden characters stop us crossing into a
    neighbouring declaration or an inline function body, while still allowing
    multi-line `static\\nReturnType foo();` declarations.
    """
    if not wanted_methods:
        return set()
    static: set[str] = set()
    for meth in wanted_methods:
        pattern = rf"\bstatic\b[^;{{}}]*\b{re.escape(meth)}\s*\("
        if re.search(pattern, class_body, re.DOTALL):
            static.add(meth)
    return static


def filter_to_relevant_methods(
    class_body: str,
    wanted_methods: set[str],
    static_methods: set[str],
) -> str:
    """
    Within a class body, mark lines that contain a wanted method name so the LLM
    knows which ones to focus on.  We don't strip other lines — the full context
    helps the LLM understand signatures — but we prepend a comment on matched lines.

    Static methods get an extra `*** STATIC ***` tag so the LLM can carry that
    fact through into the generated `Receiver:` line.

    Returns the class body with annotations added.
    """
    if not wanted_methods:
        return class_body

    lines = class_body.splitlines()
    annotated = []
    for line in lines:
        matched = [m for m in wanted_methods if m in line]
        if not matched:
            annotated.append(line)
            continue
        is_static = any(m in static_methods for m in matched)
        tag = "// *** USED IN COOT — STATIC ***" if is_static else "// *** USED IN COOT ***"
        annotated.append(f"{line}  {tag}")
    return "\n".join(annotated)


# ---------------------------------------------------------------------------
# Phase 4: Prompts & LLM calls
# ---------------------------------------------------------------------------

METHOD_DOC_SYSTEM = """\
You are an expert C++ software engineer writing concise API reference \
documentation for the MMDB2 structural biology library. \
You will be given a C++ class declaration from an MMDB2 header file. \
Lines annotated with '// *** USED IN COOT ***' are the methods you MUST \
document. Lines annotated with '// *** USED IN COOT — STATIC ***' are the \
ones declared `static` in the header — they are called via `ClassName::method(...)` \
and have no receiver. All other annotated methods are instance methods and \
must be called via an object: `obj->method(...)` or `obj.method(...)`. \
For unannotated methods you may write a one-line stub or omit them. \
Document behaviour from what the method name and signature imply — \
these are well-established structural biology data-structure methods.\
"""

METHOD_DOC_USER_TMPL = """\
Class: `{class_name}`
Methods that MUST be documented (from AST analysis of the Coot codebase):
{method_list}

Static methods (called as `{class_name}::method(...)` — no receiver):
{static_list}

All other methods in the list above are INSTANCE methods and must be called \
via an object: `obj->method(...)` or `obj.method(...)`. Never write \
`{class_name}::method(...)` for an instance method.

Here is the C++ class declaration. Lines marked '// *** USED IN COOT ***' are \
instance methods you must document. Lines marked '// *** USED IN COOT — STATIC ***' \
are static methods you must document:

```cpp
{class_body}
```

Produce documentation in this exact compact format — no markdown tables, \
no repeated field headers. Each method block looks like this:

#### `ClassName::methodName(param_type param) -> return_type`
One sentence describing what the method does.
Receiver: instance — call as `obj->methodName(...)` (or `obj.methodName(...)`).
Params: `param` — what it represents. Omit this line if there are no parameters.
Returns: what is returned and under what conditions. Omit if void.
Side effects: any state mutated or I/O performed. Omit if none.
Errors: conditions under which it fails or returns null. Omit if none.

For STATIC methods, the Receiver line MUST instead read exactly:
Receiver: static — call as `{class_name}::methodName(...)` (no receiver).

Rules:
- Use `ClassName::methodName` in the heading (always qualify with the class name).
- The `Receiver:` line is MANDATORY on every block — never omit it. It is the \
single most important field, because the model consuming this doc has been \
calling instance methods statically.
- Keep each field to one line. Be direct — no filler words.
- Omit Params/Returns/Side effects/Errors when they would be "none" or "void" \
(but never omit Receiver).
- Repeat the block for every method marked '// *** USED IN COOT ***' or \
'// *** USED IN COOT — STATIC ***'.
- Skip all unmarked methods entirely.
- Start your response with the class heading shown below.

## `{class_name}`
One sentence describing the class's overall responsibility.\
"""

INDEX_SYSTEM = """\
You are an expert C++ software engineer. Given per-class API documentation for \
the MMDB2 library as used by the Coot codebase, produce a concise API index.\
"""

INDEX_USER_TMPL = """\
Below is per-class method documentation for the MMDB2 API as used in the Coot \
codebase (refactoring target).

Produce an **API index** with:

1. **Class index** — one line per class: `ClassName (header.h)` — description.
2. **Method quick-reference** — one line per method: \
`ClassName::method(params) -> return` — one-line summary. \
Sort by class then method name.
3. **Notable patterns** — brief prose on consistent conventions: \
ownership/lifetime, error handling, naming, output parameters, etc.

---

{docs}
"""


def document_class(
    backend: LLMBackend,
    class_name: str,
    class_body: str,
    wanted_methods: set[str],
    index: int,
    total: int,
) -> str:
    static_methods = find_static_methods(class_body, wanted_methods)
    annotated = filter_to_relevant_methods(class_body, wanted_methods, static_methods)
    method_list = "\n".join(f"  - {m}" for m in sorted(wanted_methods)) or "  (all public methods)"
    static_list = (
        "\n".join(f"  - {m}" for m in sorted(static_methods))
        if static_methods
        else "  (none — every method in the list above is an instance method)"
    )

    user_prompt = METHOD_DOC_USER_TMPL.format(
        class_name=class_name,
        method_list=method_list,
        static_list=static_list,
        class_body=annotated,
    )

    print(f"  [{index}/{total}] Documenting `{class_name}` "
          f"({len(wanted_methods)} methods, {len(static_methods)} static) ...", end=" ", flush=True)
    try:
        docs = backend.complete(METHOD_DOC_SYSTEM, user_prompt, max_tokens=8192)
        print("done")
        return docs
    except Exception as e:
        print(f"error: {e}")
        return f"*Error generating documentation for `{class_name}`: {e}*"


def build_index(backend: LLMBackend, class_docs: dict[str, str]) -> str:
    parts = [f"### `{cls}`\n\n{doc}" for cls, doc in sorted(class_docs.items())]
    user_prompt = INDEX_USER_TMPL.format(docs="\n\n---\n\n".join(parts))

    print("\nGenerating API index ...", end=" ", flush=True)
    try:
        index = backend.complete(INDEX_SYSTEM, user_prompt, max_tokens=8192)
        print("done")
        return index
    except Exception as e:
        print(f"error: {e}")
        return f"*Error generating index: {e}*"


# ---------------------------------------------------------------------------
# Phase 5: Output
# ---------------------------------------------------------------------------

def build_output(class_docs: dict[str, str], api_index: str) -> str:
    lines = [
        "# MMDB2 API Reference (AST-driven)\n",
        "---\n",
        "## API Index\n",
        api_index,
        "\n---\n",
        "## Per-Class Method Documentation\n",
    ]
    for cls, doc in sorted(class_docs.items()):
        lines.append(doc)
        lines.append("\n---\n")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="AST-driven MMDB2 API documentation generator.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--report", "-r",
        required=True,
        metavar="FILE",
        help="Path to mmdb_usage_report.json (produced by the AST tool)",
    )
    parser.add_argument(
        "--mmdb-include", "-I",
        required=True,
        metavar="DIR",
        help="Path to the MMDB2 include directory containing *.h files",
    )
    parser.add_argument(
        "--output", "-o",
        default="mmdb_api_docs.md",
        help="Output markdown file (default: mmdb_api_docs.md)",
    )
    parser.add_argument(
        "--cache", "-c",
        default=None,
        metavar="FILE",
        help="JSON file to cache/resume per-class docs",
    )

    provider_group = parser.add_argument_group("provider")
    provider_group.add_argument(
        "--provider", "-p",
        choices=["anthropic", "ollama"],
        default="anthropic",
    )
    provider_group.add_argument(
        "--model", "-m",
        default=None,
        help=(
            f"Model name. "
            f"Anthropic default: {DEFAULT_ANTHROPIC_MODEL}  "
            f"Ollama default: {DEFAULT_OLLAMA_MODEL}"
        ),
    )
    provider_group.add_argument(
        "--ollama-url",
        default=DEFAULT_OLLAMA_URL,
        metavar="URL",
    )
    provider_group.add_argument(
        "--no-think",
        action="store_true",
        help=(
            "Disable model thinking/reasoning. Default is on, because the doc "
            "generator's job (classify static vs instance, infer semantics from "
            "signatures) benefits from reasoning. Use this flag to fall back to "
            "plain completions for speed or for models that don't support "
            "thinking."
        ),
    )
    parser.add_argument(
        "--skip-classes",
        nargs="*",
        default=[],
        metavar="CLASS",
        help="Class names to skip (e.g. math::Graph internal helpers)",
    )
    return parser.parse_args()


def build_backend(args: argparse.Namespace) -> LLMBackend:
    think = not args.no_think
    if args.provider == "ollama":
        return OllamaBackend(model=args.model or DEFAULT_OLLAMA_MODEL,
                             base_url=args.ollama_url,
                             think=think)
    model   = args.model or DEFAULT_ANTHROPIC_MODEL
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print("Error: ANTHROPIC_API_KEY not set.", file=sys.stderr)
        sys.exit(1)
    return AnthropicBackend(model=model, api_key=api_key, think=think)


def main() -> None:
    args    = parse_args()
    backend = build_backend(args)
    print(f"Backend      : {backend.label}")

    # --- Phase 1: Collect used symbols ---
    print(f"\nReading AST report: {args.report}")
    used = collect_used_symbols(args.report)
    print(f"Found {len(used)} unique MMDB classes referenced in the codebase")

    # Filter skipped classes
    skip = set(args.skip_classes)
    used = {cls: meths for cls, meths in used.items() if cls not in skip}

    # --- Phase 2: Map classes to headers ---
    print(f"Scanning MMDB2 headers: {args.mmdb_include}")
    class_to_header = build_class_to_header_map(args.mmdb_include)

    # Report what we can and can't find
    found_classes    = {cls for cls in used if cls in class_to_header}
    missing_classes  = {cls for cls in used if cls not in class_to_header}
    print(f"  {len(found_classes)} classes found in headers")
    if missing_classes:
        print(f"  {len(missing_classes)} classes not found in headers (will skip):")
        for cls in sorted(missing_classes):
            methods = used[cls]
            print(f"    {cls}  ({len(methods)} methods: {', '.join(sorted(methods)[:4])}{'...' if len(methods)>4 else ''})")

    # --- Load cache ---
    class_docs: dict[str, str] = {}
    if args.cache and Path(args.cache).exists():
        with open(args.cache) as f:
            class_docs = json.load(f)
        print(f"\nLoaded {len(class_docs)} cached class docs from {args.cache}")

    # --- Phase 3 + 4: Extract headers & document ---
    classes_to_process = sorted(found_classes)
    total = len(classes_to_process)
    print(f"\nDocumenting {total} classes...\n")

    for i, class_name in enumerate(classes_to_process, start=1):
        if class_name in class_docs:
            print(f"  [{i}/{total}] Skipping (cached): `{class_name}`")
            continue

        header_path = class_to_header[class_name]
        header_text = header_path.read_text(errors="replace")
        class_body  = extract_class_body(header_text, class_name)

        if not class_body:
            print(f"  [{i}/{total}] WARNING: Could not extract `{class_name}` body from {header_path.name}")
            class_docs[class_name] = f"*Could not extract class body for `{class_name}`*"
            continue

        wanted_methods = used[class_name]
        docs = document_class(backend, class_name, class_body, wanted_methods, i, total)
        class_docs[class_name] = docs

        if args.cache:
            with open(args.cache, "w") as f:
                json.dump(class_docs, f, indent=2)

        if args.provider == "anthropic":
            time.sleep(RATE_LIMIT_DELAY)

    # --- Phase 5: Index + output ---
    api_index   = build_index(backend, class_docs)
    output_path = Path(args.output)
    output_path.write_text(build_output(class_docs, api_index))
    print(f"\nAPI docs written to: {output_path}")

    # Summary
    total_methods = sum(len(v) for v in used.values() if isinstance(v, set))
    print(f"Documented {len(class_docs)} classes, {total_methods} method references")


if __name__ == "__main__":
    main()
