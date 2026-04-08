#!/usr/bin/env python3
"""
Generate per-method API documentation for a C++ codebase using an LLM.

For each C++ source/header file pair, the LLM lists every class and documents
each method: signature, parameters, return value, behaviour, side effects, and
exceptions. Then an index page is produced summarising the full public API.

Usage:
    # Anthropic Claude (default)
    python document_methods.py <root_dir> [--output api_docs.md]

    # Local Ollama model
    python document_methods.py <root_dir> --provider ollama --model codellama:34b

    # Ollama with custom host
    python document_methods.py <root_dir> --provider ollama --model llama3.1 \\
        --ollama-url http://myserver:11434
"""

import argparse
import json
import os
import sys
import time
from abc import ABC, abstractmethod
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

DEFAULT_ANTHROPIC_MODEL = "claude-opus-4-6"
DEFAULT_OLLAMA_MODEL    = "llama3.1"
DEFAULT_OLLAMA_URL      = "http://localhost:11434"
MAX_FILE_CHARS          = 80_000
RATE_LIMIT_DELAY        = 1.0


# ---------------------------------------------------------------------------
# LLM backend abstraction  (identical to summarize_codebase.py)
# ---------------------------------------------------------------------------

class LLMBackend(ABC):
    @abstractmethod
    def complete(self, system: str, user: str, max_tokens: int) -> str: ...

    @property
    @abstractmethod
    def label(self) -> str: ...


class AnthropicBackend(LLMBackend):
    def __init__(self, model: str, api_key: str) -> None:
        import anthropic
        self._client    = anthropic.Anthropic(api_key=api_key)
        self._model     = model
        self._anthropic = anthropic

    @property
    def label(self) -> str:
        return f"Claude ({self._model})"

    def complete(self, system: str, user: str, max_tokens: int) -> str:
        while True:
            try:
                with self._client.messages.stream(
                    model=self._model,
                    max_tokens=max_tokens,
                    thinking={"type": "adaptive"},
                    system=system,
                    messages=[{"role": "user", "content": user}],
                ) as stream:
                    final = stream.get_final_message()
                for block in reversed(final.content):
                    if block.type == "text":
                        return block.text
                return ""
            except self._anthropic.RateLimitError:
                print("rate limited — waiting 30s ...", end=" ", flush=True)
                time.sleep(30)


class OllamaBackend(LLMBackend):
    def __init__(self, model: str, base_url: str = DEFAULT_OLLAMA_URL) -> None:
        try:
            import httpx
        except ImportError:
            print(
                "Error: httpx is required for Ollama support.\n"
                "Install it with:  pip install httpx",
                file=sys.stderr,
            )
            sys.exit(1)

        self._httpx    = httpx
        self._model    = model
        self._base_url = base_url.rstrip("/")
        self._check_connection()

    def _check_connection(self) -> None:
        try:
            r = self._httpx.get(f"{self._base_url}/api/tags", timeout=5)
            r.raise_for_status()
        except Exception as e:
            print(
                f"Error: Cannot reach Ollama at {self._base_url}\n"
                f"  {e}\n"
                "Make sure Ollama is running (ollama serve).",
                file=sys.stderr,
            )
            sys.exit(1)

    @property
    def label(self) -> str:
        return f"Ollama ({self._model} @ {self._base_url})"

    def complete(self, system: str, user: str, max_tokens: int) -> str:
        payload = {
            "model": self._model,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user",   "content": user},
            ],
            "stream": True,
            "options": {"num_predict": max_tokens},
        }

        chunks: list[str] = []
        with self._httpx.stream(
            "POST",
            f"{self._base_url}/api/chat",
            json=payload,
            timeout=600,
        ) as response:
            response.raise_for_status()
            for line in response.iter_lines():
                if not line:
                    continue
                try:
                    data = json.loads(line)
                except json.JSONDecodeError:
                    continue
                chunk = data.get("message", {}).get("content", "")
                if chunk:
                    chunks.append(chunk)
                if data.get("done"):
                    break

        return "".join(chunks)


# ---------------------------------------------------------------------------
# File discovery & pairing  (identical to summarize_codebase.py)
# ---------------------------------------------------------------------------

SOURCE_EXTS = {".cpp", ".cc", ".cxx", ".c"}
HEADER_EXTS = {".h", ".hpp", ".hxx", ".hh"}


def find_cpp_files(root: Path) -> dict[str, dict]:
    pairs: dict[str, dict] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        ext = path.suffix.lower()
        stem_key = str(path.relative_to(root).with_suffix(""))
        if ext in SOURCE_EXTS:
            pairs.setdefault(stem_key, {})["source"] = path
        elif ext in HEADER_EXTS:
            pairs.setdefault(stem_key, {})["header"] = path
    return pairs


def read_file(path: Path) -> str:
    try:
        text = path.read_text(errors="replace")
    except OSError as e:
        return f"[Could not read {path}: {e}]"

    if len(text) > MAX_FILE_CHARS:
        half = MAX_FILE_CHARS // 2
        text = (
            text[:half]
            + f"\n\n... [truncated {len(text) - MAX_FILE_CHARS} chars] ...\n\n"
            + text[-half:]
        )
    return text


# ---------------------------------------------------------------------------
# Prompts
# ---------------------------------------------------------------------------

METHOD_DOC_SYSTEM = """\
You are an expert C++ software engineer writing precise API reference \
documentation. Your task is to read C++ source and header files and produce \
thorough, developer-facing documentation for every method in every class. \
Document behaviour from the *implementation*, not just the declaration — \
include what the method actually does internally where that is non-obvious.\
"""

METHOD_DOC_USER_TMPL = """\
Read the C++ file(s) below and produce complete method-level documentation.

For **every class and struct** found, output a section using this exact format:

---
## `ClassName`

> One-sentence description of the class's overall responsibility.

### Methods

#### `return_type methodName(param_type param_name, ...)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | One-sentence description of what this method does. |
| **Parameters** | `param_name` — type and what it represents. Repeat for each parameter. Use "none" if there are no parameters. |
| **Returns**  | What is returned and under what conditions. Use "void" if nothing is returned. |
| **Behaviour** | Step-by-step description of the internal logic (e.g. validates X, iterates Y, calls Z). |
| **Side effects** | Any state mutated, files written, signals emitted, etc. Use "none" if side-effect free. |
| **Errors / exceptions** | Conditions under which the method throws, returns an error code, or silently fails. Use "none" if not applicable. |

Repeat the method block for every method (public, protected, and private). \
Include constructors and destructors. \
If a method is only declared and not defined in these files, note that in \
**Behaviour** as "defined elsewhere". \
Do not skip any method, even trivial getters/setters.

After all classes, add a section for free functions in the same method-block \
format under a `## Free Functions` heading (omit the heading if there are none).

Use markdown throughout.

---

{file_blocks}
"""

INDEX_SYSTEM = """\
You are an expert C++ software engineer. You have per-file method documentation \
for an entire C++ codebase. Produce a concise, navigable API index.\
"""

INDEX_USER_TMPL = """\
Below is per-file method documentation for the C++ project rooted at `{root}`.

Please produce an **API index** with the following sections:

1. **Class index** – A table with columns: Class name | File | One-line \
description. Sort alphabetically by class name.

2. **Method quick-reference** – A flat table with columns: Class | Method \
signature | One-line summary. Sort by class then method name. This table is \
intended as a quick lookup for developers who know what they are looking for.

3. **Notable patterns** – A brief prose section calling out any consistent \
conventions found across the API (naming, ownership/lifetime, error handling \
strategy, use of output parameters vs return values, etc.).

Use markdown throughout.

---

{docs}
"""


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def document_file_pair(
    backend: LLMBackend,
    stem_key: str,
    pair: dict,
    index: int,
    total: int,
) -> tuple[str, str]:
    """Return (stem_key, method_docs_markdown)."""
    file_blocks_parts = []

    if "header" in pair:
        content = read_file(pair["header"])
        file_blocks_parts.append(
            f"### Header: `{pair['header'].name}`\n\n```cpp\n{content}\n```"
        )
    if "source" in pair:
        content = read_file(pair["source"])
        file_blocks_parts.append(
            f"### Source: `{pair['source'].name}`\n\n```cpp\n{content}\n```"
        )

    file_blocks = "\n\n".join(file_blocks_parts)
    user_prompt = METHOD_DOC_USER_TMPL.format(file_blocks=file_blocks)

    print(f"  [{index}/{total}] Documenting: {stem_key} ...", end=" ", flush=True)
    try:
        docs = backend.complete(METHOD_DOC_SYSTEM, user_prompt, max_tokens=8192)
        print("done")
        return stem_key, docs
    except Exception as e:
        print(f"error: {e}")
        return stem_key, f"*Error generating documentation: {e}*"


def build_index(
    backend: LLMBackend,
    root: Path,
    file_docs: dict[str, str],
) -> str:
    parts = [f"## `{k}`\n\n{v}" for k, v in sorted(file_docs.items())]
    all_docs = "\n\n---\n\n".join(parts)
    user_prompt = INDEX_USER_TMPL.format(root=root, docs=all_docs)

    print("\nGenerating API index ...", end=" ", flush=True)
    try:
        index = backend.complete(INDEX_SYSTEM, user_prompt, max_tokens=8192)
        print("done")
        return index
    except Exception as e:
        print(f"error: {e}")
        return f"*Error generating index: {e}*"


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def build_output(root: Path, file_docs: dict[str, str], api_index: str) -> str:
    lines = [
        f"# API Reference: `{root}`\n",
        "---\n",
        "## API Index\n",
        api_index,
        "\n---\n",
        "## Per-File Method Documentation\n",
    ]
    for stem_key, docs in sorted(file_docs.items()):
        lines.append(f"### `{stem_key}`\n")
        lines.append(docs)
        lines.append("\n---\n")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate per-method API docs for a C++ codebase.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("root_dir", help="Root directory of the C++ project")
    parser.add_argument(
        "--output", "-o",
        default="api_docs.md",
        help="Output markdown file (default: api_docs.md)",
    )

    provider_group = parser.add_argument_group("provider")
    provider_group.add_argument(
        "--provider", "-p",
        choices=["anthropic", "ollama"],
        default="anthropic",
        help="LLM provider (default: anthropic)",
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
        help=f"Ollama base URL (default: {DEFAULT_OLLAMA_URL})",
    )

    parser.add_argument(
        "--cache", "-c",
        default=None,
        metavar="FILE",
        help="JSON file to cache/resume per-file docs",
    )
    parser.add_argument(
        "--skip-patterns",
        nargs="*",
        default=[],
        metavar="PATTERN",
        help="Glob patterns to skip (e.g. '*/test/*' '*/third_party/*')",
    )
    return parser.parse_args()


def should_skip(stem_key: str, skip_patterns: list[str]) -> bool:
    from fnmatch import fnmatch
    return any(fnmatch(stem_key, pat) for pat in skip_patterns)


def build_backend(args: argparse.Namespace) -> LLMBackend:
    if args.provider == "ollama":
        model = args.model or DEFAULT_OLLAMA_MODEL
        return OllamaBackend(model=model, base_url=args.ollama_url)

    model   = args.model or DEFAULT_ANTHROPIC_MODEL
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print(
            "Error: ANTHROPIC_API_KEY environment variable not set.\n"
            "Set it or use --provider ollama for local inference.",
            file=sys.stderr,
        )
        sys.exit(1)
    return AnthropicBackend(model=model, api_key=api_key)


def main() -> None:
    args    = parse_args()
    root    = Path(args.root_dir).resolve()

    if not root.is_dir():
        print(f"Error: {root} is not a directory", file=sys.stderr)
        sys.exit(1)

    backend = build_backend(args)
    print(f"Backend : {backend.label}")

    print(f"Scanning {root} for C++ files ...")
    pairs = find_cpp_files(root)

    if args.skip_patterns:
        pairs = {k: v for k, v in pairs.items() if not should_skip(k, args.skip_patterns)}

    total = len(pairs)
    if total == 0:
        print("No C++ files found.")
        sys.exit(0)

    print(f"Found {total} file pair(s)/group(s).\n")

    file_docs: dict[str, str] = {}
    if args.cache and Path(args.cache).exists():
        with open(args.cache) as f:
            file_docs = json.load(f)
        print(f"Loaded {len(file_docs)} cached docs from {args.cache}\n")

    for i, (stem_key, pair) in enumerate(sorted(pairs.items()), start=1):
        if stem_key in file_docs:
            print(f"  [{i}/{total}] Skipping (cached): {stem_key}")
            continue

        _, docs = document_file_pair(backend, stem_key, pair, i, total)
        file_docs[stem_key] = docs

        if args.cache:
            with open(args.cache, "w") as f:
                json.dump(file_docs, f, indent=2)

        if args.provider == "anthropic":
            time.sleep(RATE_LIMIT_DELAY)

    api_index = build_index(backend, root, file_docs)

    output_path = Path(args.output)
    output_path.write_text(build_output(root, file_docs, api_index))
    print(f"\nAPI docs written to: {output_path}")


if __name__ == "__main__":
    main()
