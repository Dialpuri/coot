"""Report loading and global state."""
import json
import sys
from pathlib import Path


_meta: dict = {}
_files: list[dict] = []
_file_map: dict[str, dict] = {}
_coot_root: str = ""


def load_report(path: str) -> tuple[dict, list[dict], dict[str, dict], str]:
    """Load the AST report and return (meta, files, file_map, coot_root)."""
    p = Path(path)
    if not p.exists():
        print(f"ERROR: report not found at {path}", file=sys.stderr)
        sys.exit(1)
    print(f"Loading report from {path} ...")
    report = json.loads(p.read_text())
    meta     = report["meta"]
    files    = report["files"]
    file_map = {entry["rel_path"]: entry for entry in files}
    coot_root = meta.get("root", str(p.parent))
    print(f"Loaded {len(files)} files  |  {meta.get('total_mmdb_refs', 0):,} refs  |  root: {coot_root}")
    return meta, files, file_map, coot_root


def load_and_update(path: str) -> None:
    """Load the report and update module-level globals."""
    global _meta, _files, _file_map, _coot_root
    _meta, _files, _file_map, _coot_root = load_report(path)
