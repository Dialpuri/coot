import argparse
import json
import os
import sys
from pathlib import Path

import httpx
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

# ── Config ──────────────────────────────────────────────────────────────────
_parser = argparse.ArgumentParser(description="MMDB Refactor Manager backend")
_parser.add_argument(
    "--report",
    default="/Users/dialpuri/lmb/coot/mmdb_usage_report.json",
    help="Path to mmdb_usage_report.json produced by the AST tool",
)
_args, _ = _parser.parse_known_args()

REPORT_PATH = _args.report
PROGRESS_PATH = Path(__file__).parent / "progress.json"
OLLAMA_URL = "http://localhost:11434/api/generate"

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


# ── Helpers ──────────────────────────────────────────────────────────────────
def load_progress() -> dict:
    if PROGRESS_PATH.exists():
        with open(PROGRESS_PATH, "r") as f:
            return json.load(f)
    return {}


def save_progress(data: dict) -> None:
    with open(PROGRESS_PATH, "w") as f:
        json.dump(data, f, indent=2)


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

    return StreamingResponse(stream_ollama(), media_type="text/plain")


@app.get("/api/progress")
def get_progress():
    return load_progress()


@app.post("/api/progress")
def update_progress(update: ProgressUpdate):
    data = load_progress()
    data[update.key] = update.status
    save_progress(data)
    return {"ok": True, "key": update.key, "status": update.status}
