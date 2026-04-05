#!/usr/bin/env python3
"""
visualise_ast_report.py

Reads mmdb_usage_report.json (produced by the AST tool) and writes a
self-contained interactive HTML dashboard.

Usage:
    python3 mmdb-recon/visualise_ast_report.py [--input FILE] [--output FILE]
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict


# ---------------------------------------------------------------------------
# Data transforms  (all run at generation time, no source code embedded)
# ---------------------------------------------------------------------------

def build_stats(meta, files):
    total_funcs = sum(len(f["functions"]) for f in files)
    return {
        "total_files_scanned": meta.get("total_files_scanned", 0),
        "files_with_mmdb":     meta.get("files_with_mmdb_refs", 0),
        "total_refs":          meta.get("total_mmdb_refs", 0),
        "total_funcs":         total_funcs,
    }


def build_file_rows(files):
    rows = []
    for f in files:
        rows.append({
            "path":   f["rel_path"],
            "dir":    f["rel_path"].split("/")[0],
            "refs":   f["total_mmdb_refs"],
            "funcs":  len(f["functions"]),
        })
    return sorted(rows, key=lambda r: -r["refs"])


def build_function_index(files):
    """
    Returns a dict:  rel_path -> list of function dicts.
    Source code is NOT embedded — the detail pane shows mmdb refs only.
    """
    index = {}
    for f in files:
        funcs = []
        for fn in f["functions"]:
            funcs.append({
                "name":      fn["name"],
                "line":      fn["line"],
                "end_line":  fn.get("end_line", fn["line"]),
                "refs":      fn["mmdb_ref_count"],
                "symbols":   fn["mmdb_symbols"],
                "mmdb_refs": fn.get("refs", []),
            })
        index[f["rel_path"]] = funcs
    return index


def build_symbol_freq(files):
    counts = defaultdict(int)
    for f in files:
        for fn in f["functions"]:
            for sym in fn["mmdb_symbols"]:
                counts[sym] += fn["mmdb_ref_count"]
    return sorted(
        [{"symbol": s, "refs": r} for s, r in counts.items()],
        key=lambda x: -x["refs"]
    )[:40]


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>MMDB AST Report</title>
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:system-ui,-apple-system,sans-serif;background:#f0f2f5;color:#1a1a2e;height:100vh;display:flex;flex-direction:column}

/* ── Header ── */
.header{background:#1a1a2e;color:#fff;padding:12px 20px;flex-shrink:0}
.header h1{font-size:1rem;font-weight:600;letter-spacing:.4px;margin-bottom:8px}
.stats{display:flex;gap:16px;flex-wrap:wrap}
.stat{background:rgba(255,255,255,.1);border-radius:5px;padding:6px 14px}
.stat .val{font-size:1.3rem;font-weight:700;color:#64dfdf}
.stat .lbl{font-size:.68rem;opacity:.65;text-transform:uppercase;letter-spacing:.5px}

/* ── Layout ── */
.body{display:grid;grid-template-columns:340px 1fr;flex:1;overflow:hidden}

/* ── Left pane: file list ── */
.pane-left{background:#fff;border-right:1px solid #e0e0e0;display:flex;flex-direction:column;overflow:hidden}
.pane-left .toolbar{padding:8px 10px;border-bottom:1px solid #eee;display:flex;gap:6px;flex-shrink:0}
.pane-left .toolbar input{flex:1;border:1px solid #ddd;border-radius:4px;padding:5px 8px;font-size:.8rem;outline:none}
.pane-left .toolbar input:focus{border-color:#4361ee}
.pane-left .toolbar select{border:1px solid #ddd;border-radius:4px;padding:5px 6px;font-size:.78rem;outline:none;max-width:110px}
.file-list{flex:1;overflow-y:auto}
.file-item{padding:7px 12px;border-bottom:1px solid #f2f2f2;cursor:pointer;display:flex;align-items:baseline;gap:6px;transition:background .1s}
.file-item:hover{background:#f5f7ff}
.file-item.active{background:#eef2ff;border-left:3px solid #4361ee}
.file-item .name{font-size:.78rem;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:#333}
.file-item .dir-tag{font-size:.68rem;color:#888;background:#f0f0f0;border-radius:3px;padding:1px 5px;flex-shrink:0}
.file-item .refs-badge{font-size:.72rem;font-weight:700;color:#4361ee;flex-shrink:0;min-width:30px;text-align:right}
.file-item .fns-badge{font-size:.68rem;color:#999;flex-shrink:0;min-width:30px;text-align:right}
.file-count{font-size:.72rem;color:#888;padding:5px 12px;border-bottom:1px solid #eee;background:#fafafa;flex-shrink:0}

/* ── Right pane ── */
.pane-right{display:flex;flex-direction:column;overflow:hidden}

/* ── Function list ── */
.fn-pane{flex:1;overflow-y:auto;padding:12px 16px;min-height:0}
.fn-pane h2{font-size:.8rem;font-weight:700;text-transform:uppercase;letter-spacing:.5px;color:#888;margin-bottom:10px}
.fn-toolbar{display:flex;gap:6px;margin-bottom:10px;flex-wrap:wrap}
.fn-toolbar input{border:1px solid #ddd;border-radius:4px;padding:5px 8px;font-size:.8rem;outline:none;width:220px}
.fn-toolbar input:focus{border-color:#4361ee}
.fn-table{width:100%;border-collapse:collapse;font-size:.8rem}
.fn-table th{text-align:left;padding:6px 8px;border-bottom:2px solid #e0e0e0;font-size:.72rem;font-weight:700;color:#555;text-transform:uppercase;letter-spacing:.3px;cursor:pointer;white-space:nowrap}
.fn-table th:hover{color:#4361ee}
.fn-table td{padding:6px 8px;border-bottom:1px solid #f0f0f0;vertical-align:top}
.fn-table tr:hover td{background:#f9faff}
.fn-table tr.fn-active td{background:#eef2ff}
.fn-name{font-weight:600;color:#1a1a2e;cursor:pointer}
.fn-name:hover{color:#4361ee;text-decoration:underline}
.fn-loc{color:#999;font-size:.72rem}
.sym-chips{display:flex;flex-wrap:wrap;gap:3px;margin-top:3px}
.chip{background:#e8eaf6;color:#3949ab;border-radius:3px;padding:1px 5px;font-size:.68rem}
.chip.highlight{background:#fff3cd;color:#856404}

/* ── Detail drawer (bottom) ── */
.detail-pane{border-top:2px solid #e0e0e0;background:#fff;display:flex;flex-direction:column;overflow:hidden;transition:height .2s;height:0}
.detail-pane.open{height:300px}
.detail-header{display:flex;align-items:center;padding:8px 14px;border-bottom:1px solid #eee;flex-shrink:0;gap:10px}
.detail-header .fn-title{font-size:.85rem;font-weight:700;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.detail-header .loc-tag{font-size:.72rem;color:#64dfdf;background:#1a1a2e;border-radius:4px;padding:2px 8px;flex-shrink:0}
.detail-header .close-btn{cursor:pointer;color:#aaa;font-size:1.1rem;line-height:1;padding:2px 6px;border-radius:4px}
.detail-header .close-btn:hover{background:#f0f0f0;color:#333}
.ref-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:3px}
.ref-item{padding:4px 8px;font-size:.75rem;display:flex;gap:8px;align-items:baseline;background:#f8f9ff;border-radius:3px}
.ref-item .ref-line{color:#4361ee;font-weight:700;min-width:40px;font-family:monospace;flex-shrink:0}
.ref-item .ref-sym{color:#333;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}

/* ── Symbol panel ── */
.sym-panel{border-top:1px solid #e0e0e0;background:#fff;padding:10px 16px;flex-shrink:0}
.sym-panel h2{font-size:.72rem;font-weight:700;text-transform:uppercase;letter-spacing:.4px;color:#aaa;margin-bottom:8px}
.sym-bar-list{display:flex;flex-wrap:wrap;gap:5px}

/* ── Empty state ── */
.empty{display:flex;align-items:center;justify-content:center;height:100%;color:#bbb;font-size:.9rem;flex-direction:column;gap:8px}
.empty svg{opacity:.3}

/* ── Scrollbars ── */
::-webkit-scrollbar{width:6px;height:6px}
::-webkit-scrollbar-track{background:transparent}
::-webkit-scrollbar-thumb{background:#ccc;border-radius:3px}
</style>
</head>
<body>

<div class="header">
  <h1>MMDB AST Usage Report &mdash; Refactoring Dashboard</h1>
  <div class="stats" id="stats"></div>
</div>

<div class="body">
  <!-- Left: file list -->
  <div class="pane-left">
    <div class="toolbar">
      <input id="file-search" placeholder="Search files…" oninput="filterFiles()">
      <select id="dir-filter" onchange="filterFiles()">
        <option value="">All dirs</option>
      </select>
    </div>
    <div class="file-count" id="file-count"></div>
    <div class="file-list" id="file-list"></div>
  </div>

  <!-- Right: functions + detail -->
  <div class="pane-right">
    <div class="fn-pane" id="fn-pane">
      <div class="empty">
        <svg width="40" height="40" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
          <path d="M3 7h18M3 12h18M3 17h18"/>
        </svg>
        Select a file to see its functions
      </div>
    </div>

    <!-- Detail drawer -->
    <div class="detail-pane" id="detail-pane">
      <div class="detail-header">
        <span class="fn-title" id="detail-title"></span>
        <span class="loc-tag" id="detail-loc"></span>
        <span class="close-btn" onclick="closeDetail()">✕</span>
      </div>
      <div id="refs-col" style="overflow-y:auto;padding:8px 14px;flex:1"></div>
    </div>

    <!-- Symbol freq bar -->
    <div class="sym-panel">
      <h2>Top mmdb symbols by reference count</h2>
      <div class="sym-bar-list" id="sym-bar"></div>
    </div>
  </div>
</div>

<script>
// ── Embedded data ──
const STATS    = __STATS__;
const FILES    = __FILES__;
const FN_INDEX = __FN_INDEX__;
const SYM_FREQ = __SYM_FREQ__;

// ── State ──
let activeFile = null;
let activeFn   = null;
let fnSortCol  = "refs";
let fnSortDir  = -1;

// ── Init ──
(function init() {
  // Stats bar
  const s = STATS;
  document.getElementById("stats").innerHTML = [
    { lbl: "Files scanned",    val: s.total_files_scanned.toLocaleString() },
    { lbl: "Files with MMDB",  val: s.files_with_mmdb.toLocaleString() },
    { lbl: "Total refs",       val: s.total_refs.toLocaleString() },
    { lbl: "Functions",        val: s.total_funcs.toLocaleString() },
  ].map(c => `<div class="stat"><div class="val">${c.val}</div><div class="lbl">${c.lbl}</div></div>`).join("");

  // Dir filter options
  const dirs = [...new Set(FILES.map(f => f.dir))].sort();
  const sel = document.getElementById("dir-filter");
  dirs.forEach(d => sel.add(new Option(d, d)));

  renderFileList(FILES);
  renderSymBar();
})();

// ── File list ──
function filterFiles() {
  const q   = document.getElementById("file-search").value.toLowerCase();
  const dir = document.getElementById("dir-filter").value;
  const filtered = FILES.filter(f =>
    (!q   || f.path.toLowerCase().includes(q)) &&
    (!dir || f.dir === dir)
  );
  renderFileList(filtered);
}

function renderFileList(rows) {
  document.getElementById("file-count").textContent =
    `${rows.length} file${rows.length !== 1 ? "s" : ""}`;

  document.getElementById("file-list").innerHTML = rows.map(f => `
    <div class="file-item ${activeFile === f.path ? "active" : ""}"
         onclick="selectFile('${f.path.replace(/'/g, "\\'")}')">
      <span class="name" title="${f.path}">${f.path}</span>
      <span class="dir-tag">${f.dir}</span>
      <span class="fns-badge">${f.funcs} fn</span>
      <span class="refs-badge">${f.refs}</span>
    </div>`).join("");
}

// ── Function list ──
function selectFile(path) {
  activeFile = path;
  activeFn   = null;
  closeDetail();

  // Re-render file list to update active state
  filterFiles();

  const fns = FN_INDEX[path] || [];
  const pane = document.getElementById("fn-pane");

  if (fns.length === 0) {
    pane.innerHTML = `<div class="empty">No functions found in this file</div>`;
    return;
  }

  pane.innerHTML = `
    <h2>${path}</h2>
    <div class="fn-toolbar">
      <input id="fn-search" placeholder="Filter functions…" oninput="renderFnTable()">
    </div>
    <table class="fn-table" id="fn-table">
      <thead>
        <tr>
          <th onclick="sortFns('name')">Function</th>
          <th onclick="sortFns('refs')" style="width:60px;text-align:right">Refs</th>
          <th onclick="sortFns('line')" style="width:80px">Lines</th>
        </tr>
      </thead>
      <tbody id="fn-tbody"></tbody>
    </table>`;

  renderFnTable();
}

function renderFnTable() {
  const q    = (document.getElementById("fn-search")?.value || "").toLowerCase();
  const fns  = (FN_INDEX[activeFile] || [])
    .filter(fn => !q || fn.name.toLowerCase().includes(q))
    .slice().sort((a, b) => fnSortDir * (a[fnSortCol] < b[fnSortCol] ? -1 : a[fnSortCol] > b[fnSortCol] ? 1 : 0));

  document.getElementById("fn-tbody").innerHTML = fns.map(fn => `
    <tr class="${activeFn === fn.name + fn.line ? "fn-active" : ""}"
        onclick="selectFn('${activeFile.replace(/'/g, "\\'")}', ${fn.line})">
      <td>
        <div class="fn-name">${escHtml(fn.name)}</div>
        <div class="sym-chips">${fn.symbols.slice(0, 6).map(s =>
          `<span class="chip">${escHtml(s.replace("class ",""))}</span>`
        ).join("")}${fn.symbols.length > 6 ? `<span class="chip">+${fn.symbols.length - 6}</span>` : ""}</div>
      </td>
      <td style="text-align:right;font-weight:700;color:#4361ee">${fn.refs}</td>
      <td class="fn-loc">${fn.line}–${fn.end_line}</td>
    </tr>`).join("");
}

function sortFns(col) {
  fnSortDir = (fnSortCol === col) ? -fnSortDir : -1;
  fnSortCol = col;
  renderFnTable();
}

// ── Detail drawer ──
function selectFn(filePath, line) {
  const fns = FN_INDEX[filePath] || [];
  const fn  = fns.find(f => f.line === line);
  if (!fn) return;

  activeFn = fn.name + fn.line;
  renderFnTable();

  // Header
  document.getElementById("detail-title").textContent = fn.name;
  document.getElementById("detail-loc").textContent =
    `${filePath}  :${fn.line}–${fn.end_line}  (${fn.refs} refs)`;

  // Refs grid
  const refsHtml = fn.mmdb_refs.length
    ? `<div style="font-size:.72rem;font-weight:700;text-transform:uppercase;letter-spacing:.4px;color:#aaa;margin-bottom:8px">
         MMDB References — ${fn.mmdb_refs.length} occurrences
       </div>
       <div class="ref-grid">` +
      fn.mmdb_refs.map(([l, sym]) =>
        `<div class="ref-item">
           <span class="ref-line">:${l}</span>
           <span class="ref-sym" title="${escHtml(sym)}">${escHtml(sym.replace("class ",""))}</span>
         </div>`).join("") + `</div>`
    : `<div style="color:#bbb;font-size:.8rem;padding:8px 0">No inline refs recorded</div>`;

  document.getElementById("refs-col").innerHTML = refsHtml;

  // Open drawer
  document.getElementById("detail-pane").classList.add("open");
}

function closeDetail() {
  activeFn = null;
  document.getElementById("detail-pane").classList.remove("open");
}

// ── Symbol bar ──
function renderSymBar() {
  const max = SYM_FREQ[0]?.refs || 1;
  document.getElementById("sym-bar").innerHTML = SYM_FREQ.map(s => {
    const pct = Math.round(100 * s.refs / max);
    const name = s.symbol.replace("class ", "").replace("mmdb::", "");
    return `<div title="${escHtml(s.symbol)}: ${s.refs} refs"
      style="background:linear-gradient(90deg,#c7d2fe ${pct}%,#e8eaf6 ${pct}%);
             border-radius:4px;padding:3px 8px;font-size:.72rem;color:#3949ab;cursor:default;white-space:nowrap">
      ${escHtml(name)} <strong>${s.refs}</strong>
    </div>`;
  }).join("");
}

// ── Utility ──
function escHtml(s) {
  return String(s)
    .replace(/&/g,"&amp;").replace(/</g,"&lt;")
    .replace(/>/g,"&gt;").replace(/"/g,"&quot;");
}
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Build + main
# ---------------------------------------------------------------------------

def build_html(stats, file_rows, fn_index, sym_freq):
    def j(obj):
        return json.dumps(obj, separators=(",", ":")).replace("</", "<\\/")

    html = HTML
    html = html.replace("__STATS__",    j(stats))
    html = html.replace("__FILES__",    j(file_rows))
    html = html.replace("__FN_INDEX__", j(fn_index))
    html = html.replace("__SYM_FREQ__", j(sym_freq))
    return html


def main():
    parser = argparse.ArgumentParser(description="Visualise AST mmdb report as HTML.")
    parser.add_argument("--input",  default="mmdb_usage_report.json")
    parser.add_argument("--output", default="mmdb_ast_report.html")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: {input_path} not found", file=sys.stderr)
        sys.exit(1)

    print(f"Reading {input_path} …", file=sys.stderr)
    report    = json.loads(input_path.read_text())
    meta      = report["meta"]
    files     = report["files"]
    coot_root = meta.get("root", ".")

    print("Building data …", file=sys.stderr)
    stats    = build_stats(meta, files)
    rows     = build_file_rows(files)
    fn_index = build_function_index(files)
    sym_freq = build_symbol_freq(files)

    print("Generating HTML …", file=sys.stderr)
    html = build_html(stats, rows, fn_index, sym_freq)

    out = Path(args.output)
    out.write_text(html)
    print(f"Written to {out}  ({out.stat().st_size // 1024} KB)", file=sys.stderr)


if __name__ == "__main__":
    main()
