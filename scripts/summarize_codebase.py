#!/usr/bin/env python3
"""
Summarize a C++ codebase using an LLM (Claude or local Ollama).

For each C++ source/header file pair, generates a class-level summary.
Then produces an overall project summary from all file summaries.

Usage:
    # Anthropic Claude (default)
    python summarize_codebase.py <root_dir> [--output summary.md]

    # Local Ollama model
    python summarize_codebase.py <root_dir> --provider ollama --model codellama:34b

    # Ollama with custom host
    python summarize_codebase.py <root_dir> --provider ollama --model llama3.1 \\
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
MAX_FILE_CHARS          = 80_000   # ~20K tokens; truncate larger files
RATE_LIMIT_DELAY        = 1.0      # seconds between Anthropic API calls


# ---------------------------------------------------------------------------
# LLM backend abstraction
# ---------------------------------------------------------------------------

class LLMBackend(ABC):
    """Minimal interface: given a system prompt + user message, return text."""

    @abstractmethod
    def complete(self, system: str, user: str, max_tokens: int) -> str: ...

    @property
    @abstractmethod
    def label(self) -> str: ...


class AnthropicBackend(LLMBackend):
    def __init__(self, model: str, api_key: str) -> None:
        import anthropic
        self._client = anthropic.Anthropic(api_key=api_key)
        self._model  = model
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
                # Return the last text block (thinking blocks come first)
                for block in reversed(final.content):
                    if block.type == "text":
                        return block.text
                return ""
            except self._anthropic.RateLimitError:
                print("rate limited — waiting 30s ...", end=" ", flush=True)
                time.sleep(30)


class OllamaBackend(LLMBackend):
    """Calls a local Ollama instance via its native /api/chat endpoint.

    Uses httpx (already a dependency of the anthropic package) so no
    extra install is required.
    """

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
        """Send a chat request and stream the response back."""
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
            timeout=600,          # large files can take a while locally
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
# File discovery & pairing
# ---------------------------------------------------------------------------

SOURCE_EXTS = {".cpp", ".cc", ".cxx", ".c"}
HEADER_EXTS = {".h", ".hpp", ".hxx", ".hh"}


def find_cpp_files(root: Path) -> dict[str, dict]:
    """Return pairs keyed by stem (relative path without extension)."""
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

FILE_SUMMARY_SYSTEM = """\
You are an expert C++ software engineer performing a codebase audit. \
Your task is to read C++ source and header files and produce a concise, \
structured technical summary that will later be used to build an overall \
project understanding. Focus on what the code *does*, not how it is formatted.\
"""

FILE_SUMMARY_USER_TMPL = """\
Please analyse the following C++ file(s) and produce a structured summary \
covering:

1. **Purpose** – What problem does this file/class solve? What is its role \
in the broader system?
2. **Key classes / structs** – For each class or struct: name, responsibility, \
important data members, and notable methods.
3. **Key free functions** – Any significant non-member functions and what they do.
4. **Dependencies** – What external libraries, headers, or project modules does \
this file rely on?
5. **Notable design patterns or algorithms** – Any non-obvious techniques worth \
calling out (e.g. RAII, visitor pattern, custom allocator, specific algorithm).
6. **Potential issues or TODOs** – Anything that looks incomplete, fragile, or \
flagged with TODO/FIXME comments.

Keep the summary concise but technically precise. Use markdown.

---

{file_blocks}
"""

OVERALL_SUMMARY_SYSTEM = """\
You are an expert C++ software engineer. You have been given per-file technical \
summaries of an entire C++ codebase. Your task is to synthesise these into a \
single, coherent project-level overview.\
"""

OVERALL_SUMMARY_USER_TMPL = """\
Below are per-file summaries of every C++ source and header file in the project \
rooted at `{root}`.

Please write a comprehensive **project overview** that covers:

1. **Project purpose** – What does this codebase do at a high level?
2. **Architecture** – How is the code organised? What are the major subsystems \
or modules and how do they interact?
3. **Key abstractions** – The most important classes, interfaces, or data \
structures that a new developer must understand first.
4. **Data flow** – How does data enter the system, get processed, and exit?
5. **External dependencies** – Third-party libraries or frameworks relied upon.
6. **Build / entry points** – Any apparent main() functions, library entry \
points, or public APIs.
7. **Suggested reading order** – For a developer new to this codebase, which \
files should they read first and why?

Use markdown with clear headings.

---

{summaries}
"""


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def summarise_file_pair(
    backend: LLMBackend,
    stem_key: str,
    pair: dict,
    index: int,
    total: int,
) -> tuple[str, str]:
    """Return (stem_key, summary_markdown)."""
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
    user_prompt = FILE_SUMMARY_USER_TMPL.format(file_blocks=file_blocks)

    print(f"  [{index}/{total}] Summarising: {stem_key} ...", end=" ", flush=True)
    try:
        summary = backend.complete(FILE_SUMMARY_SYSTEM, user_prompt, max_tokens=4096)
        print("done")
        return stem_key, summary
    except Exception as e:
        print(f"error: {e}")
        return stem_key, f"*Error generating summary: {e}*"


def summarise_project(
    backend: LLMBackend,
    root: Path,
    file_summaries: dict[str, str],
) -> str:
    parts = [f"## `{k}`\n\n{v}" for k, v in sorted(file_summaries.items())]
    all_summaries = "\n\n---\n\n".join(parts)
    user_prompt = OVERALL_SUMMARY_USER_TMPL.format(root=root, summaries=all_summaries)

    print("\nGenerating overall project summary ...", end=" ", flush=True)
    try:
        summary = backend.complete(OVERALL_SUMMARY_SYSTEM, user_prompt, max_tokens=8192)
        print("done")
        return summary
    except Exception as e:
        print(f"error: {e}")
        return f"*Error generating project summary: {e}*"


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def build_output(root: Path, file_summaries: dict[str, str], project_summary: str) -> str:
    lines = [
        f"# Codebase Summary: `{root}`\n",
        "---\n",
        "## Project Overview\n",
        project_summary,
        "\n---\n",
        "## Per-File Summaries\n",
    ]
    for stem_key, summary in sorted(file_summaries.items()):
        lines.append(f"### `{stem_key}`\n")
        lines.append(summary)
        lines.append("\n---\n")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarise a C++ codebase using an LLM.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("root_dir", help="Root directory of the C++ project")
    parser.add_argument(
        "--output", "-o",
        default="codebase_summary.md",
        help="Output markdown file (default: codebase_summary.md)",
    )

    # Provider selection
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

    # Misc
    parser.add_argument(
        "--cache", "-c",
        default=None,
        metavar="FILE",
        help="JSON file to cache/resume per-file summaries",
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

    # Anthropic
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
    args   = parse_args()
    root   = Path(args.root_dir).resolve()

    if not root.is_dir():
        print(f"Error: {root} is not a directory", file=sys.stderr)
        sys.exit(1)

    backend = build_backend(args)
    print(f"Backend : {backend.label}")

    # --- Discover files ---
    print(f"Scanning {root} for C++ files ...")
    pairs = find_cpp_files(root)

    if args.skip_patterns:
        pairs = {k: v for k, v in pairs.items() if not should_skip(k, args.skip_patterns)}

    total = len(pairs)
    if total == 0:
        print("No C++ files found.")
        sys.exit(0)

    print(f"Found {total} file pair(s)/group(s).\n")

    # --- Load cache ---
    file_summaries: dict[str, str] = {}
    if args.cache and Path(args.cache).exists():
        with open(args.cache) as f:
            file_summaries = json.load(f)
        print(f"Loaded {len(file_summaries)} cached summaries from {args.cache}\n")

    # --- Summarise each file pair ---
    for i, (stem_key, pair) in enumerate(sorted(pairs.items()), start=1):
        if stem_key in file_summaries:
            print(f"  [{i}/{total}] Skipping (cached): {stem_key}")
            continue

        _, summary = summarise_file_pair(backend, stem_key, pair, i, total)
        file_summaries[stem_key] = summary

        if args.cache:
            with open(args.cache, "w") as f:
                json.dump(file_summaries, f, indent=2)

        # Anthropic rate-limit courtesy delay; Ollama doesn't need it
        if args.provider == "anthropic":
            time.sleep(RATE_LIMIT_DELAY)

    # --- Overall project summary ---
    project_summary = summarise_project(backend, root, file_summaries)

    # --- Write output ---
    output_path = Path(args.output)
    output_path.write_text(build_output(root, file_summaries, project_summary))
    print(f"\nSummary written to: {output_path}")


if __name__ == "__main__":
    main()
