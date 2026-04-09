"""Endpoints: /api/stats, /api/files, /api/file/{path}, /api/reload, /api/source"""
import os

from fastapi import APIRouter, HTTPException, Query

import report
from config import REPORT_PATH

router = APIRouter()


@router.get("/api/stats")
def get_stats():
    return {
        "total_files_scanned":  report._meta.get("total_files_scanned", 0),
        "files_with_mmdb_refs": report._meta.get("files_with_mmdb_refs", 0),
        "total_mmdb_refs":      report._meta.get("total_mmdb_refs", 0),
    }


@router.get("/api/files")
def get_files(
    search: str = Query(""),
    sort: str = Query("refs"),
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
):
    items = report._files

    if search:
        search_lower = search.lower()
        items = [f for f in items if search_lower in f["rel_path"].lower()]

    if sort == "refs":
        items = sorted(items, key=lambda f: f.get("total_mmdb_refs", 0), reverse=True)
    elif sort == "path":
        items = sorted(items, key=lambda f: f["rel_path"])
    elif sort == "functions":
        items = sorted(items, key=lambda f: len(f.get("functions", [])), reverse=True)

    total = len(items)
    start = (page - 1) * page_size
    page_items = items[start: start + page_size]

    return {
        "items": [
            {
                "rel_path":        f["rel_path"],
                "total_mmdb_refs": f.get("total_mmdb_refs", 0),
                "function_count":  len(f.get("functions", [])),
                "includes_mmdb":   f.get("includes_mmdb", False),
            }
            for f in page_items
        ],
        "total":     total,
        "page":      page,
        "page_size": page_size,
    }


@router.get("/api/file/{path:path}")
def get_file(path: str):
    entry = report._file_map.get(path)
    if entry is None:
        raise HTTPException(status_code=404, detail=f"File not found in report: {path}")
    return entry


@router.post("/api/reload")
def reload_report():
    """Reload the report from disk without restarting the server."""
    report.load_and_update(REPORT_PATH)
    return {
        "ok":              True,
        "files":           len(report._files),
        "total_mmdb_refs": report._meta.get("total_mmdb_refs", 0),
        "root":            report._coot_root,
    }


@router.get("/api/source")
def get_source(
    path: str = Query(...),
    start: int = Query(1, ge=1),
    end: int = Query(0),
):
    full_path = os.path.join(report._coot_root, path)
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
        "end":   s + len(selected),
    }
