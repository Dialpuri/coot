#!/bin/bash
# Scans all .cc/.cpp files in the coot source tree (excluding dot-folders and
# build directories) for mmdb namespace usage using the AST tool.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="$SCRIPT_DIR/build/find_mmdb_functions"
COOT_ROOT="/Users/dialpuri/lmb/coot"
COMPILE_DB="$COOT_ROOT/cmake-build-debug"
OUTPUT="$COOT_ROOT/mmdb_usage_report.json"
BATCH_DIR="$(mktemp -d)"
BATCH_SIZE=50

if [ ! -f "$TOOL" ]; then
  echo "ERROR: tool not found at $TOOL — run 'make' in $SCRIPT_DIR/build first"
  exit 1
fi

echo "Collecting source files..."
FILE_LIST=$(find "$COOT_ROOT" -type f \( -name "*.cc" -o -name "*.cpp" \) \
  | grep -v '/\.' \
  | grep -v '/cmake-build' \
  | sort)

TOTAL=$(echo "$FILE_LIST" | wc -l | tr -d ' ')
echo "Found $TOTAL source files"

# Split into batches and run the tool on each
BATCH_NUM=0
BATCH_FILES=()
BATCH_OUTPUTS=()

run_batch() {
  if [ ${#BATCH_FILES[@]} -eq 0 ]; then return; fi
  BATCH_OUT="$BATCH_DIR/batch_${BATCH_NUM}.json"
  BATCH_OUTPUTS+=("$BATCH_OUT")
  echo "  Running batch $BATCH_NUM (${#BATCH_FILES[@]} files)..."
  "$TOOL" \
    -p "$COMPILE_DB" \
    --output="$BATCH_OUT" \
    "${BATCH_FILES[@]}" 2>/dev/null || true
  BATCH_NUM=$((BATCH_NUM + 1))
  BATCH_FILES=()
}

COUNT=0
while IFS= read -r f; do
  BATCH_FILES+=("$f")
  COUNT=$((COUNT + 1))
  if [ ${#BATCH_FILES[@]} -ge $BATCH_SIZE ]; then
    run_batch
  fi
  if [ $((COUNT % 100)) -eq 0 ]; then
    echo "  Queued $COUNT / $TOTAL files..."
  fi
done <<< "$FILE_LIST"

# Run any remaining files
run_batch

echo "Ran $BATCH_NUM batches. Merging results..."

# Merge all batch JSON files into a single report
python3 - "$OUTPUT" "${BATCH_OUTPUTS[@]}" <<'PYEOF'
import json, sys
from pathlib import Path
from collections import defaultdict

out_path = sys.argv[1]
batch_files = sys.argv[2:]

all_files_map = {}   # rel_path -> file entry
total_files_scanned = 0
coot_root = None

for bf in batch_files:
    p = Path(bf)
    if not p.exists():
        continue
    data = json.loads(p.read_text())
    if coot_root is None:
        coot_root = data.get("meta", {}).get("root", ".")
    total_files_scanned += data.get("meta", {}).get("total_files_scanned", 0)

    for entry in data.get("files", []):
        rel = entry["rel_path"]
        if rel not in all_files_map:
            all_files_map[rel] = entry
        else:
            # Merge functions if the same file appeared in multiple batches
            existing_names = {f["name"] for f in all_files_map[rel]["functions"]}
            for fn in entry.get("functions", []):
                if fn["name"] not in existing_names:
                    all_files_map[rel]["functions"].append(fn)
                    all_files_map[rel]["total_mmdb_refs"] += fn["mmdb_ref_count"]

merged_files = sorted(all_files_map.values(), key=lambda f: -f["total_mmdb_refs"])
total_mmdb_refs = sum(f["total_mmdb_refs"] for f in merged_files)

report = {
    "meta": {
        "root": coot_root or ".",
        "total_files_scanned": total_files_scanned,
        "files_with_mmdb_refs": len(merged_files),
        "total_mmdb_refs": total_mmdb_refs,
    },
    "files": merged_files,
}

Path(out_path).write_text(json.dumps(report, indent=2))
print(f"Merged {len(batch_files)} batches -> {len(merged_files)} files with mmdb refs, {total_mmdb_refs} total refs")
print(f"Report written to: {out_path}")
PYEOF

# Cleanup temp batch files
rm -rf "$BATCH_DIR"

echo "Done. Regenerate the HTML dashboard with:"
echo "  python3 $COOT_ROOT/mmdb-recon/visualise_mmdb_usage.py --input $OUTPUT"
