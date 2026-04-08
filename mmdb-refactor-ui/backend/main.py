import argparse
import asyncio
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import httpx
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

import doc_extractor

# ── Config ──────────────────────────────────────────────────────────────────
_parser = argparse.ArgumentParser(description="MMDB Refactor Manager backend")
_parser.add_argument(
    "--report",
    default="../mmdb-recon/mmdb_usage_report.json",
    help="Path to mmdb_usage_report.json produced by the AST tool",
)
_parser.add_argument(
    "--cxx", default="c++",
    help="C++ compiler to use for test compilation",
)
_parser.add_argument(
    "--cxx-flags", default="",
    help="Extra flags to pass to the compiler (e.g. -I/extra/include)",
)
_parser.add_argument(
    "--mmdb-docs",
    default="../mmdb-recon/llm/mmdb_methods.md",
    help="Path to the mmdb_methods.md file produced by document_methods.py",
)
_args, _ = _parser.parse_known_args()

REPORT_PATH = _args.report
PROGRESS_PATH = Path(__file__).parent / "progress.json"
TESTS_PATH = Path(__file__).parent / "tests.json"
OLLAMA_URL = "http://localhost:11434/api/generate"

# ── MMDB API docs (loaded once at startup) ───────────────────────────────────
_mmdb_docs_markdown = doc_extractor.load_docs(_args.mmdb_docs)
if _mmdb_docs_markdown:
    print(f"Loaded MMDB docs from {_args.mmdb_docs} ({len(_mmdb_docs_markdown):,} chars)")
else:
    print(f"WARNING: MMDB docs not found at {_args.mmdb_docs} — prompts will lack API context")


# ── GTest detection ───────────────────────────────────────────────────────────
def _detect_gtest_flags() -> str:
    """Try pkg-config first, then fall back to common install prefixes."""
    try:
        flags = subprocess.check_output(
            ["pkg-config", "--cflags", "--libs", "gtest", "gtest_main"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        if flags:
            return flags
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass
    for prefix in ["/opt/homebrew", "/usr/local", "/usr"]:
        if Path(f"{prefix}/include/gtest/gtest.h").exists():
            return f"-I{prefix}/include -L{prefix}/lib -lgtest -lgtest_main"
    return "-lgtest -lgtest_main"


_gtest_flags = _detect_gtest_flags()
print(f"GTest flags: {_gtest_flags}")

# ── Strategy context injected into every refactor prompt ────────────────────
STRATEGY_CONTEXT = """You are an expert C++ developer helping to migrate code from the MMDB2 structural biology library to the Gemmi library.

## Type Mapping: MMDB → Gemmi

### Atom Properties
| MMDB | Gemmi |
|------|-------|
| atom->x, atom->y, atom->z | atom.pos.x, atom.pos.y, atom.pos.z |
| atom->occupancy | atom.occ |
| atom->tempFactor | atom.b_iso |
| atom->element (const char*) | atom.element.name() |
| atom->name | atom.name |
| atom->serNum | atom.serial |
| atom->Distance(other) | atom.pos.distance(other.pos) |

### Residue Information
| MMDB | Gemmi |
|------|-------|
| res->GetResName() | res.name |
| res->GetSeqNum() | res.seqid.num |
| res->GetInsCode() | res.seqid.icode |
| res->GetNumberOfAtoms() | res.atoms.size() |

### Pointer Types
| MMDB | Gemmi |
|------|-------|
| mmdb::Atom* / mmdb::PAtom | gemmi::Atom& or gemmi::Atom* |
| mmdb::PPAtom (array) | std::vector<gemmi::Atom*> |
| mmdb::Residue* / mmdb::PResidue | gemmi::Residue& or gemmi::Residue* |
| mmdb::PPResidue | std::vector<gemmi::Residue*> |

### Enums / Selection Flags
| MMDB | Gemmi approach |
|------|----------------|
| STYPE_RESIDUE, STYPE_ATOM | Use container iteration |
| ANY_RES | No range check needed |
| SKEY_NEW | Not needed |
| mmdb::realtype | double |
| mmdb::mat44 | gemmi::Transform |

### Math
| MMDB | Gemmi |
|------|-------|
| mmdb::math::Alignment | Custom or Gemmi alignment |
| mmdb::mat44 | gemmi::Transform |
| mmdb::rvector | gemmi::Vec3 |
| atom->Transform(m) | atom.pos = transform.apply(atom.pos) |

### File Loading
MMDB Pattern:
  mmdb::Manager *mol = new mmdb::Manager();
  mol->ReadStructure("structure.pdb");
Gemmi Equivalent:
  gemmi::Structure st = gemmi::read_structure("structure.pdb");
  gemmi::Model *model = &st.models[0];

### Hierarchy Navigation
MMDB iterates with index-based GetModel/GetChain/GetResidue/GetAtom loops.
Gemmi uses range-for over st.models, model.chains, chain.residues, residue.atoms.

### Atom Selection/Queries
MMDB uses NewSelection/Select/GetSelIndex/DeleteSelection with handle integers.
Gemmi uses direct iteration with conditional filtering on chain.name, res.seqid, etc.

Replace all mmdb:: types, pointer arithmetic, and selection handles with idiomatic Gemmi code.
Prefer references over raw pointers where possible. Remove manual memory management."""

# ── Report loader ────────────────────────────────────────────────────────────
def load_report(path: str) -> tuple[dict, list[dict], dict[str, dict], str]:
    """Load the AST report and return (meta, files, file_map, coot_root)."""
    p = Path(path)
    if not p.exists():
        print(f"ERROR: report not found at {path}", file=sys.stderr)
        sys.exit(1)
    print(f"Loading report from {path} ...")
    report = json.loads(p.read_text())
    meta      = report["meta"]
    files     = report["files"]
    file_map  = {entry["rel_path"]: entry for entry in files}
    coot_root = meta.get("root", str(p.parent))
    print(f"Loaded {len(files)} files  |  {meta.get('total_mmdb_refs', 0):,} refs  |  root: {coot_root}")
    return meta, files, file_map, coot_root

_meta, _files, _file_map, _coot_root = load_report(REPORT_PATH)

# ── App ──────────────────────────────────────────────────────────────────────
app = FastAPI(title="MMDB Refactor Manager")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── Models ──────────────────────────────────────────────────────────────────
class RefactorRequest(BaseModel):
    function_name: str
    code: str
    mmdb_symbols: list[str]
    model: str = "codellama"
    additional_instructions: str = ""


class ProgressUpdate(BaseModel):
    key: str
    status: str


class TestRecord(BaseModel):
    key: str
    mmdb_test: str = ""
    gemmi_test: str = ""
    notes: str = ""
    status: str = "draft"  # draft | reviewed | done


class GenerateTestRequest(BaseModel):
    function_name: str
    source_code: str
    mmdb_symbols: list[str]
    target: str = "both"   # "mmdb" | "gemmi" | "both"
    model: str = "gemma4"
    additional_instructions: str = ""


class WriteTestRequest(BaseModel):
    rel_source_path: str   # e.g. "coot-utils/coot-coord-utils.cc"
    fn_name: str           # e.g. "coot::util::delete_residue_references_in_header_info"
    fn_line: int
    mmdb_test: str = ""
    gemmi_test: str = ""


class CompileRunRequest(BaseModel):
    rel_source_path: str
    fn_name: str
    fn_line: int
    variant: str           # "mmdb" | "gemmi"
    test_code: str         # current textarea content (written before compile)


class GitCommitRequest(BaseModel):
    rel_source_path: str
    fn_name: str
    variant: str           # "mmdb" | "gemmi" | "both"
    commit_message: str


class GenerateAllRequest(BaseModel):
    model: str = "gemma4"
    skip_existing: bool = True
    additional_instructions: str = ""


class ValidateFixRequest(BaseModel):
    rel_source_path: str
    fn_name: str
    fn_line: int
    variant: str        # "mmdb" | "gemmi"
    test_code: str
    model: str = "gemma4"


# ── Helpers ──────────────────────────────────────────────────────────────────
def load_progress() -> dict:
    if PROGRESS_PATH.exists():
        with open(PROGRESS_PATH, "r") as f:
            return json.load(f)
    return {}


def save_progress(data: dict) -> None:
    with open(PROGRESS_PATH, "w") as f:
        json.dump(data, f, indent=2)


def load_tests() -> dict:
    if TESTS_PATH.exists():
        with open(TESTS_PATH, "r") as f:
            return json.load(f)
    return {}


def save_tests(data: dict) -> None:
    with open(TESTS_PATH, "w") as f:
        json.dump(data, f, indent=2)


def sanitize_fn_name(name: str) -> str:
    """coot::util::foo_bar  ->  coot_util_foo_bar"""
    return re.sub(r'[^a-zA-Z0-9]', '_', name).strip('_')


def get_test_file_path(rel_source_path: str, fn_name: str, variant: str) -> Path:
    """Return the absolute path where a test .cc file should be written."""
    source_dir = Path(rel_source_path).parent          # e.g.  coot-utils
    tests_dir  = Path(_coot_root) / source_dir / "tests"
    safe       = sanitize_fn_name(fn_name)
    return tests_dir / f"{safe}_{variant}_test.cc"


_VARIANT_INCLUDES: dict[str, list[str]] = {
    "mmdb":  ["#include <mmdb2/mmdb_manager.h>"],
    "gemmi": ["#include <gemmi/structure.hpp>", "#include <gemmi/model.hpp>"],
}

MAX_TEST_RETRIES = 5

# ── Compile/run constants (shared with batch retry) ───────────────────────────
_COOT_API_DIR  = "/Users/dialpuri/lmb/build-coot-and-deps/"
_COOT_API_NAME = "cootapi"
_MMDB_API_DIR  = "/opt/homebrew/Cellar/mmdb2/2.0.22/lib"
_MMDB_API_NAME = "mmdb2"


def wrap_test_content(content: str, variant: str, fn_name: str) -> str:
    """Ensure the .cc file has a gtest header, variant includes, and a main()."""
    content = content.strip()
    header: list[str] = [f"// Auto-generated {variant.upper()} test for {fn_name}"]
    if "#include <gtest/gtest.h>" not in content:
        header.append("#include <gtest/gtest.h>")
    for inc in _VARIANT_INCLUDES.get(variant, []):
        if inc not in content:
            header.append(inc)
    needs_main = "RUN_ALL_TESTS()" not in content and "int main(" not in content
    parts = ["\n".join(header), "", content]
    if needs_main:
        parts += ["", "int main(int argc, char **argv) {",
                  "  ::testing::InitGoogleTest(&argc, argv);",
                  "  return RUN_ALL_TESTS();", "}"]
    return "\n".join(parts) + "\n"


# ── Endpoints ────────────────────────────────────────────────────────────────

@app.get("/api/stats")
def get_stats():
    return {
        "total_files_scanned": _meta.get("total_files_scanned", 0),
        "files_with_mmdb_refs": _meta.get("files_with_mmdb_refs", 0),
        "total_mmdb_refs": _meta.get("total_mmdb_refs", 0),
    }


@app.get("/api/files")
def get_files(
    search: str = Query(""),
    sort: str = Query("refs"),
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
):
    items = _files

    # Filter
    if search:
        search_lower = search.lower()
        items = [f for f in items if search_lower in f["rel_path"].lower()]

    # Sort
    if sort == "refs":
        items = sorted(items, key=lambda f: f.get("total_mmdb_refs", 0), reverse=True)
    elif sort == "path":
        items = sorted(items, key=lambda f: f["rel_path"])
    elif sort == "functions":
        items = sorted(items, key=lambda f: len(f.get("functions", [])), reverse=True)

    total = len(items)
    start = (page - 1) * page_size
    end = start + page_size
    page_items = items[start:end]

    return {
        "items": [
            {
                "rel_path": f["rel_path"],
                "total_mmdb_refs": f.get("total_mmdb_refs", 0),
                "function_count": len(f.get("functions", [])),
                "includes_mmdb": f.get("includes_mmdb", False),
            }
            for f in page_items
        ],
        "total": total,
        "page": page,
        "page_size": page_size,
    }


@app.get("/api/file/{path:path}")
def get_file(path: str):
    entry = _file_map.get(path)
    if entry is None:
        raise HTTPException(status_code=404, detail=f"File not found in report: {path}")
    return entry


@app.post("/api/reload")
def reload_report():
    """Reload the report from disk without restarting the server."""
    global _meta, _files, _file_map, _coot_root
    _meta, _files, _file_map, _coot_root = load_report(REPORT_PATH)
    return {
        "ok": True,
        "files": len(_files),
        "total_mmdb_refs": _meta.get("total_mmdb_refs", 0),
        "root": _coot_root,
    }


@app.get("/api/source")
def get_source(
    path: str = Query(...),
    start: int = Query(1, ge=1),
    end: int = Query(0),
):
    full_path = os.path.join(_coot_root, path)
    try:
        with open(full_path, "r", errors="replace") as f:
            all_lines = f.readlines()
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Source file not found: {full_path}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading file: {e}")

    total = len(all_lines)
    s = max(0, start - 1)
    e = min(total, end if end > 0 else total)
    selected = all_lines[s:e]

    return {
        "lines": "".join(selected),
        "start": s + 1,
        "end": s + len(selected),
    }


@app.post("/api/refactor")
async def refactor(req: RefactorRequest):
    symbols_list = ", ".join(req.mmdb_symbols) if req.mmdb_symbols else "none listed"

    prompt = (
        f"Refactor the following C++ function from MMDB2 to Gemmi.\n\n"
        f"Function name: {req.function_name}\n"
        f"MMDB symbols used: {symbols_list}\n"
    )
    if req.additional_instructions.strip():
        prompt += f"\nAdditional instructions: {req.additional_instructions.strip()}\n"

    prompt += (
        f"\n## Original MMDB code:\n```cpp\n{req.code}\n```\n\n"
        "## Task:\n"
        "Rewrite the function above using Gemmi equivalents. "
        "Output only the refactored C++ code in a single ```cpp``` block. "
        "Do not include explanation outside the code block unless it is a brief inline comment."
    )

    payload = {
        "model": req.model,
        "system": STRATEGY_CONTEXT,
        "prompt": prompt,
        "stream": True,
    }

    async def stream_ollama():
        async with httpx.AsyncClient(timeout=300.0) as client:
            async with client.stream("POST", OLLAMA_URL, json=payload) as response:
                if response.status_code != 200:
                    body = await response.aread()
                    yield f"Error from Ollama ({response.status_code}): {body.decode()}"
                    return
                async for line in response.aiter_lines():
                    if line.strip():
                        try:
                            data = json.loads(line)
                            chunk = data.get("response", "")
                            if chunk:
                                yield chunk
                        except json.JSONDecodeError:
                            pass

    return StreamingResponse(stream_ollama(), media_type="text/plain")  # /api/refactor


@app.get("/api/progress")
def get_progress():
    return load_progress()


@app.post("/api/progress")
def update_progress(update: ProgressUpdate):
    data = load_progress()
    data[update.key] = update.status
    save_progress(data)
    return {"ok": True, "key": update.key, "status": update.status}


# ── Test endpoints ────────────────────────────────────────────────────────────

@app.get("/api/tests")
def get_tests():
    return load_tests()


@app.post("/api/tests")
def save_test(record: TestRecord):
    data = load_tests()
    data[record.key] = {
        "mmdb_test": record.mmdb_test,
        "gemmi_test": record.gemmi_test,
        "notes": record.notes,
        "status": record.status,
    }
    save_tests(data)
    return {"ok": True, "key": record.key}


@app.get("/api/tests/load")
def load_test_for_function(
    rel_source_path: str = Query(...),
    fn_name: str = Query(...),
    fn_line: int = Query(...),
):
    """Return the test record for a specific function, falling back to .cc files on disk."""
    key = f"{rel_source_path}::{fn_name}:{fn_line}"
    record = load_tests().get(key, {})

    mmdb_content = record.get("mmdb_test", "")
    gemmi_content = record.get("gemmi_test", "")

    # Fallback: read directly from the .cc file if JSON entry is absent/empty
    if not mmdb_content:
        mmdb_path = get_test_file_path(rel_source_path, fn_name, "mmdb")
        if mmdb_path.exists():
            mmdb_content = mmdb_path.read_text(errors="replace")

    if not gemmi_content:
        gemmi_path = get_test_file_path(rel_source_path, fn_name, "gemmi")
        if gemmi_path.exists():
            gemmi_content = gemmi_path.read_text(errors="replace")

    return {
        "mmdb_test": mmdb_content,
        "gemmi_test": gemmi_content,
        "notes": record.get("notes", ""),
        "status": record.get("status", "draft"),
    }


def _build_api_context(symbols: list[str]) -> str:
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


@app.post("/api/generate-test")
async def generate_test(req: GenerateTestRequest):
    symbols_list = ", ".join(req.mmdb_symbols) if req.mmdb_symbols else "none"
    api_context  = _build_api_context(req.mmdb_symbols)

    if req.target == "mmdb":
        prompt = (
            f"Write a Google Test (GTest) unit test for the following C++ function "
            f"using MMDB2 types.\n\n"
            f"Function: `{req.function_name}`\n"
            f"MMDB symbols used: {symbols_list}\n"
        )
        if req.additional_instructions.strip():
            prompt += f"Additional instructions: {req.additional_instructions.strip()}\n"
        prompt += (
            f"{api_context}\n\n"
            f"## Source code:\n```cpp\n{req.source_code}\n```\n\n"
            "Write a focused GTest unit test that:\n"
            "- Sets up minimal MMDB2 structures needed to call the function\n"
            "- Calls the function with realistic inputs\n"
            "- Asserts the expected outputs/side-effects\n"
            "- Uses TEST() or TEST_F() macros\n\n"
            "Output ONLY the test code inside a single ```cpp``` block."
        )
        system = STRATEGY_CONTEXT

    elif req.target == "gemmi":
        prompt = (
            f"Write a Google Test (GTest) unit test for the following C++ function "
            f"AFTER it has been refactored from MMDB2 to Gemmi.\n\n"
            f"Function: `{req.function_name}`\n"
            f"Original MMDB symbols: {symbols_list}\n"
        )
        if req.additional_instructions.strip():
            prompt += f"Additional instructions: {req.additional_instructions.strip()}\n"
        prompt += (
            f"{api_context}\n\n"
            f"## Original MMDB source (to be refactored):\n```cpp\n{req.source_code}\n```\n\n"
            "Write a GTest unit test for the Gemmi-refactored version that:\n"
            "- Sets up gemmi::Structure / gemmi::Model / gemmi::Chain / gemmi::Residue / gemmi::Atom\n"
            "- Tests the SAME logical behaviour as the MMDB version would\n"
            "- Uses TEST() or TEST_F() macros\n\n"
            "Output ONLY the test code inside a single ```cpp``` block."
        )
        system = STRATEGY_CONTEXT

    else:  # both
        prompt = (
            f"Write TWO Google Test (GTest) unit tests for the following C++ function:\n"
            f"1. An MMDB2 version (tests the current implementation)\n"
            f"2. A Gemmi version (tests the future refactored implementation)\n\n"
            f"Function: `{req.function_name}`\n"
            f"MMDB symbols used: {symbols_list}\n"
        )
        if req.additional_instructions.strip():
            prompt += f"Additional instructions: {req.additional_instructions.strip()}\n"
        prompt += (
            f"{api_context}\n\n"
            f"## Source code:\n```cpp\n{req.source_code}\n```\n\n"
            "Both tests must verify the SAME logical behaviour with equivalent assertions, "
            "but use each library's own types and APIs for setup.\n\n"
            "Output EXACTLY this format — no other text:\n\n"
            "=== MMDB TEST ===\n"
            "```cpp\n"
            "// MMDB2 version\n"
            "TEST(FunctionNameTests, MMDB_Behaviour) {\n"
            "  // ...\n"
            "}\n"
            "```\n\n"
            "=== GEMMI TEST ===\n"
            "```cpp\n"
            "// Gemmi version\n"
            "TEST(FunctionNameTests, Gemmi_Behaviour) {\n"
            "  // ...\n"
            "}\n"
            "```"
        )
        system = STRATEGY_CONTEXT

    payload = {
        "model": req.model,
        "system": system,
        "prompt": prompt,
        "stream": True,
    }

    async def stream_ollama():
        async with httpx.AsyncClient(timeout=300.0) as client:
            async with client.stream("POST", OLLAMA_URL, json=payload) as response:
                if response.status_code != 200:
                    body = await response.aread()
                    yield f"Error from Ollama ({response.status_code}): {body.decode()}"
                    return
                async for line in response.aiter_lines():
                    if line.strip():
                        try:
                            data = json.loads(line)
                            chunk = data.get("response", "")
                            if chunk:
                                yield chunk
                        except json.JSONDecodeError:
                            pass

    return StreamingResponse(stream_ollama(), media_type="text/plain")  # /api/generate-test


# ── Compile/run helpers ───────────────────────────────────────────────────────

def _make_compile_cmd(test_file: Path, output_bin: Path) -> str:
    return (
        f'{_args.cxx} -std=c++17 "{test_file}" -o "{output_bin}" '
        f'{_gtest_flags} -I"{_coot_root}" -pthread '
        f'-Wl,-rpath,{_COOT_API_DIR} -L {_COOT_API_DIR} -L {_MMDB_API_DIR} '
        f'-l{_COOT_API_NAME} -l {_MMDB_API_NAME}'
    )


async def _compile_test(test_file: Path) -> tuple[bool, str]:
    """Compile a test file; return (success, output)."""
    bin_dir = test_file.parent / "bin"
    bin_dir.mkdir(exist_ok=True)
    output_bin = bin_dir / test_file.stem
    cmd = _make_compile_cmd(test_file, output_bin)
    if _args.cxx_flags:
        cmd += f" {_args.cxx_flags}"
    proc = await asyncio.create_subprocess_shell(
        cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT,
        cwd=str(_coot_root),
    )
    stdout, _ = await proc.communicate()
    return proc.returncode == 0, stdout.decode(errors="replace")


async def _run_test(test_file: Path) -> tuple[bool, str]:
    """Run a compiled test binary; return (success, output)."""
    output_bin = test_file.parent / "bin" / test_file.stem
    if not output_bin.exists():
        return False, "binary not found"
    run_proc = await asyncio.create_subprocess_exec(
        str(output_bin), stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT,
    )
    stdout, _ = await run_proc.communicate()
    return run_proc.returncode == 0, stdout.decode(errors="replace")


def _build_fix_prompt(fn_name: str, variant: str, current_code: str, error: str, attempt: int) -> str:
    inc_hint = (
        "Required headers:\n#include <gtest/gtest.h>\n#include <mmdb2/mmdb_manager.h>"
        if variant == "mmdb" else
        "Required headers:\n#include <gtest/gtest.h>\n#include <gemmi/structure.hpp>\n#include <gemmi/model.hpp>"
    )
    return (
        f"A C++ GTest failed to compile or run (attempt {attempt} of {MAX_TEST_RETRIES}). Fix it.\n\n"
        f"Function: `{fn_name}`  |  Variant: {variant.upper()}\n\n"
        f"{inc_hint}\n\n"
        f"## Current code:\n```cpp\n{current_code}\n```\n\n"
        f"## Error:\n```\n{error[:2000]}\n```\n\n"
        "Output ONLY the corrected test in a single ```cpp``` block. "
        "Include all necessary headers at the top."
    )


# ── Test file write / compile / run endpoints ─────────────────────────────────

@app.get("/api/tests/file-paths")
def get_test_file_paths(rel_source_path: str = Query(...), fn_name: str = Query(...)):
    """Return where MMDB and Gemmi test files would be written."""
    mmdb_path = get_test_file_path(rel_source_path, fn_name, "mmdb")
    gemmi_path = get_test_file_path(rel_source_path, fn_name, "gemmi")
    return {
        "mmdb": str(mmdb_path.relative_to(_coot_root)),
        "gemmi": str(gemmi_path.relative_to(_coot_root)),
        "mmdb_exists": mmdb_path.exists(),
        "gemmi_exists": gemmi_path.exists(),
    }


@app.post("/api/tests/write")
def write_test_files(req: WriteTestRequest):
    """Write MMDB and/or Gemmi test .cc files to disk."""
    written: dict[str, str] = {}
    for variant, content in [("mmdb", req.mmdb_test), ("gemmi", req.gemmi_test)]:
        if not content.strip():
            continue
        path = get_test_file_path(req.rel_source_path, req.fn_name, variant)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(wrap_test_content(content, variant, req.fn_name))
        written[variant] = str(path.relative_to(_coot_root))
    return {"ok": True, "written": written}


@app.post("/api/tests/compile-run")
async def compile_and_run(req: CompileRunRequest):
    """Write the test file, compile it, run it, stream all output."""
    test_file = get_test_file_path(req.rel_source_path, req.fn_name, req.variant)
    test_file.parent.mkdir(parents=True, exist_ok=True)
    test_file.write_text(wrap_test_content(req.test_code, req.variant, req.fn_name))

    rel_path = test_file.relative_to(_coot_root)
    bin_dir = test_file.parent / "bin"
    bin_dir.mkdir(exist_ok=True)
    output_bin = bin_dir / test_file.stem

    compile_cmd = _make_compile_cmd(test_file, output_bin)
    if _args.cxx_flags:
        compile_cmd += f" {_args.cxx_flags}"

    async def stream_output():
        yield f"[WRITE]   {rel_path}\n"
        yield f"[COMPILE] {compile_cmd}\n\n"
        yield f"[WORKING_DIR] {str(_coot_root)}\n\n"

        proc = await asyncio.create_subprocess_shell(
            compile_cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            cwd=str(_coot_root),
        )
        assert proc.stdout is not None
        async for raw in proc.stdout:
            yield raw.decode(errors="replace")
        await proc.wait()

        if proc.returncode != 0:
            yield f"\n[COMPILE FAILED]  exit code {proc.returncode}\n"
            return

        yield f"\n[COMPILE OK]\n[RUN]     {output_bin.relative_to(_coot_root)}\n\n"

        run_proc = await asyncio.create_subprocess_exec(
            str(output_bin),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
        )
        assert run_proc.stdout is not None
        async for raw in run_proc.stdout:
            yield raw.decode(errors="replace")
        await run_proc.wait()

        code = run_proc.returncode
        yield f"\n{'[PASSED]' if code == 0 else f'[FAILED]  exit code {code}'}\n"

    return StreamingResponse(stream_output(), media_type="text/plain")  # compile-run


@app.post("/api/tests/git-commit")
async def git_commit_test(req: GitCommitRequest):
    """Git-add and commit test files. The caller must confirm before invoking."""
    variants = ["mmdb", "gemmi"] if req.variant == "both" else [req.variant]
    files_to_add = []
    for v in variants:
        p = get_test_file_path(req.rel_source_path, req.fn_name, v)
        if p.exists():
            files_to_add.append(str(p))

    async def stream_output():
        if not files_to_add:
            yield "[GIT ERROR] No existing test files found to commit — write them first\n"
            return

        rel_paths = [str(Path(f).relative_to(_coot_root)) for f in files_to_add]
        yield f"[GIT ADD]  {', '.join(rel_paths)}\n"

        add_proc = await asyncio.create_subprocess_shell(
            "git add " + " ".join(f'"{f}"' for f in files_to_add),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            cwd=str(_coot_root),
        )
        assert add_proc.stdout is not None
        async for raw in add_proc.stdout:
            yield raw.decode(errors="replace")
        await add_proc.wait()

        if add_proc.returncode != 0:
            yield f"[GIT ADD FAILED]  exit code {add_proc.returncode}\n"
            return

        # Escape any double-quotes in the message
        safe_msg = req.commit_message.replace('"', "'")
        commit_proc = await asyncio.create_subprocess_shell(
            f'git commit -m "{safe_msg}"',
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            cwd=str(_coot_root),
        )
        assert commit_proc.stdout is not None
        async for raw in commit_proc.stdout:
            yield raw.decode(errors="replace")
        await commit_proc.wait()

        if commit_proc.returncode != 0:
            yield f"\n[GIT COMMIT FAILED]  exit code {commit_proc.returncode}\n"
        else:
            yield "\n[COMMITTED]\n"

    return StreamingResponse(stream_output(), media_type="text/plain")  # git-commit


# ── Batch test generation helpers ─────────────────────────────────────────────

def _parse_both_sections_py(raw: str) -> tuple[str, str]:
    mmdb_marker = "=== MMDB TEST ==="
    gemmi_marker = "=== GEMMI TEST ==="
    mi, gi = raw.find(mmdb_marker), raw.find(gemmi_marker)
    if mi == -1 and gi == -1:
        return raw, ""
    if mi != -1 and gi == -1:
        return raw[mi + len(mmdb_marker):].strip(), ""
    if mi == -1:
        return "", raw[gi + len(gemmi_marker):].strip()
    return raw[mi + len(mmdb_marker):gi].strip(), raw[gi + len(gemmi_marker):].strip()


def _strip_fences_py(code: str) -> str:
    m = re.match(r'^```(?:cpp)?\s*\n([\s\S]*?)```\s*$', code.strip())
    return m.group(1) if m else code


def _build_both_prompt(fn_name: str, source: str, symbols: list[str], extra: str) -> str:
    symbols_list = ", ".join(symbols) if symbols else "none"
    prompt = (
        f"Write TWO Google Test (GTest) unit tests for the following C++ function:\n"
        f"1. An MMDB2 version (tests the current implementation)\n"
        f"2. A Gemmi version (tests the future refactored implementation)\n\n"
        f"Function: `{fn_name}`\n"
        f"MMDB symbols used: {symbols_list}\n"
    )
    if extra.strip():
        prompt += f"Additional instructions: {extra.strip()}\n"
    prompt += (
        f"\n## Source code:\n```cpp\n{source}\n```\n\n"
        "Both tests must verify the SAME logical behaviour with equivalent assertions, "
        "but use each library's own types and APIs for setup.\n\n"
        "Required headers for each variant:\n"
        "  MMDB test:  #include <gtest/gtest.h>  #include <mmdb2/mmdb_manager.h>\n"
        "  Gemmi test: #include <gtest/gtest.h>  #include <gemmi/structure.hpp>  #include <gemmi/model.hpp>\n\n"
        "Output EXACTLY this format — no other text:\n\n"
        "=== MMDB TEST ===\n"
        "```cpp\n"
        "// MMDB2 version\n"
        "TEST(FunctionNameTests, MMDB_Behaviour) {\n"
        "  // ...\n"
        "}\n"
        "```\n\n"
        "=== GEMMI TEST ===\n"
        "```cpp\n"
        "// Gemmi version\n"
        "TEST(FunctionNameTests, Gemmi_Behaviour) {\n"
        "  // ...\n"
        "}\n"
        "```"
    )
    return prompt


async def _call_ollama(model: str, prompt: str) -> str:
    """Non-streaming Ollama call — waits for full response."""
    payload = {"model": model, "system": STRATEGY_CONTEXT, "prompt": prompt, "stream": False}
    async with httpx.AsyncClient(timeout=600.0) as client:
        resp = await client.post(OLLAMA_URL, json=payload)
        resp.raise_for_status()
        return resp.json().get("response", "")


@app.post("/api/tests/validate-fix")
async def validate_and_fix(req: ValidateFixRequest):
    """
    Compile + run a single test variant. On failure, ask the LLM to fix it and
    retry up to MAX_TEST_RETRIES times. Streams NDJSON progress events.

    Event types: start | compiling | compile_output | running | run_output |
                 fixing | fixed_code | done
    """
    async def stream():
        path = get_test_file_path(req.rel_source_path, req.fn_name, req.variant)
        path.parent.mkdir(parents=True, exist_ok=True)
        current_code = req.test_code

        yield json.dumps({"type": "start", "max": MAX_TEST_RETRIES, "variant": req.variant}) + "\n"

        for attempt in range(1, MAX_TEST_RETRIES + 1):
            path.write_text(wrap_test_content(current_code, req.variant, req.fn_name))

            yield json.dumps({"type": "compiling", "attempt": attempt, "max": MAX_TEST_RETRIES}) + "\n"
            compile_ok, compile_out = await _compile_test(path)
            yield json.dumps({"type": "compile_output", "text": compile_out, "ok": compile_ok}) + "\n"

            if compile_ok:
                yield json.dumps({"type": "running", "attempt": attempt}) + "\n"
                run_ok, run_out = await _run_test(path)
                yield json.dumps({"type": "run_output", "text": run_out, "ok": run_ok}) + "\n"
                if run_ok:
                    yield json.dumps({"type": "done", "status": "pass", "attempts": attempt, "code": current_code}) + "\n"
                    return
                error_out = f"Compile: OK\nRun failed:\n{run_out}"
            else:
                error_out = compile_out

            if attempt < MAX_TEST_RETRIES:
                yield json.dumps({"type": "fixing", "attempt": attempt}) + "\n"
                fix_prompt = _build_fix_prompt(req.fn_name, req.variant, current_code, error_out, attempt)
                fixed_raw = await _call_ollama(req.model, fix_prompt)
                current_code = _strip_fences_py(fixed_raw) or current_code
                yield json.dumps({"type": "fixed_code", "code": current_code}) + "\n"
            else:
                yield json.dumps({"type": "done", "status": "fail", "attempts": attempt, "code": current_code}) + "\n"

    return StreamingResponse(stream(), media_type="application/x-ndjson")


@app.post("/api/tests/generate-all")
async def generate_all_tests(req: GenerateAllRequest):
    """
    Iterate every function, generate both tests, compile+run each, retry up to
    MAX_TEST_RETRIES times on failure using the compile error as LLM feedback.
    Streams NDJSON: start | progress | skip | attempt | done | error | finish
    """
    async def stream_progress():
        tests = load_tests()
        work = [
            (entry["rel_path"], fn)
            for entry in _files
            for fn in entry.get("functions", [])
        ]
        total = len(work)
        yield json.dumps({"type": "start", "total": total}) + "\n"

        done = skipped = errors = 0

        for rel_path, fn in work:
            fn_name: str = fn["name"]
            fn_line: int = fn["line"]
            fn_end: int  = fn["end_line"]
            key = f"{rel_path}::{fn_name}:{fn_line}"

            if req.skip_existing:
                rec = tests.get(key, {})
                if rec.get("mmdb_test") and rec.get("gemmi_test"):
                    skipped += 1
                    yield json.dumps({
                        "type": "skip", "key": key, "fn": fn_name, "file": rel_path,
                        "done": done, "skipped": skipped, "errors": errors, "total": total,
                    }) + "\n"
                    continue

            yield json.dumps({
                "type": "progress", "key": key, "fn": fn_name, "file": rel_path,
                "done": done, "skipped": skipped, "errors": errors, "total": total,
            }) + "\n"

            try:
                # ── Load source ──────────────────────────────────────────────
                full_path = os.path.join(_coot_root, rel_path)
                with open(full_path, "r", errors="replace") as f:
                    all_lines = f.readlines()
                source = "".join(all_lines[max(0, fn_line - 1):min(len(all_lines), fn_end)])

                # ── Initial LLM generation ───────────────────────────────────
                prompt = _build_both_prompt(fn_name, source, fn.get("mmdb_symbols", []), req.additional_instructions)
                raw = await _call_ollama(req.model, prompt)
                mmdb_raw, gemmi_raw = _parse_both_sections_py(raw)
                initial = {"mmdb": _strip_fences_py(mmdb_raw), "gemmi": _strip_fences_py(gemmi_raw)}

                # ── Compile-run-retry loop per variant ───────────────────────
                variant_results: dict[str, dict] = {}
                for variant in ("mmdb", "gemmi"):
                    current_code = initial[variant]
                    if not current_code.strip():
                        variant_results[variant] = {"code": "", "status": "skip", "attempts": 0}
                        continue

                    path = get_test_file_path(rel_path, fn_name, variant)
                    path.parent.mkdir(parents=True, exist_ok=True)

                    for attempt in range(1, MAX_TEST_RETRIES + 1):
                        path.write_text(wrap_test_content(current_code, variant, fn_name))
                        compile_ok, compile_out = await _compile_test(path)

                        if compile_ok:
                            run_ok, run_out = await _run_test(path)
                            if run_ok:
                                variant_results[variant] = {"code": current_code, "status": "pass", "attempts": attempt}
                                break
                            error_out = f"Compile: OK\nRun output:\n{run_out}"
                        else:
                            error_out = compile_out

                        if attempt < MAX_TEST_RETRIES:
                            yield json.dumps({
                                "type": "attempt",
                                "key": key, "fn": fn_name, "file": rel_path,
                                "variant": variant, "attempt": attempt, "max": MAX_TEST_RETRIES,
                                "error": error_out[:600],
                                "done": done, "skipped": skipped, "errors": errors, "total": total,
                            }) + "\n"
                            fix_prompt = _build_fix_prompt(fn_name, variant, current_code, error_out, attempt)
                            fixed_raw = await _call_ollama(req.model, fix_prompt)
                            current_code = _strip_fences_py(fixed_raw) or current_code
                        else:
                            variant_results[variant] = {"code": current_code, "status": "fail", "attempts": attempt}

                # ── Persist best result to tests.json ────────────────────────
                rec = tests.get(key, {"notes": "", "status": "draft"})
                for variant, result in variant_results.items():
                    if result["code"].strip():
                        rec[f"{variant}_test"] = result["code"]
                tests[key] = rec
                save_tests(tests)

                mmdb_r = variant_results.get("mmdb",  {"status": "skip", "attempts": 0})
                gemmi_r = variant_results.get("gemmi", {"status": "skip", "attempts": 0})
                done += 1
                yield json.dumps({
                    "type": "done", "key": key, "fn": fn_name, "file": rel_path,
                    "done": done, "skipped": skipped, "errors": errors, "total": total,
                    "mmdb_status":   mmdb_r["status"],   "mmdb_attempts":  mmdb_r["attempts"],
                    "gemmi_status":  gemmi_r["status"],  "gemmi_attempts": gemmi_r["attempts"],
                }) + "\n"

            except Exception as ex:
                errors += 1
                yield json.dumps({
                    "type": "error", "key": key, "fn": fn_name, "file": rel_path,
                    "done": done, "skipped": skipped, "errors": errors, "total": total,
                    "message": str(ex),
                }) + "\n"

        yield json.dumps({
            "type": "finish", "done": done, "skipped": skipped, "errors": errors, "total": total,
        }) + "\n"

    return StreamingResponse(stream_progress(), media_type="application/x-ndjson")
