import { useState, useEffect, useRef, useCallback } from 'react'
import type { FileDetail, FunctionRecord } from '../types'
import { fetchSource, saveTest, fetchTestFilePaths, writeTestFiles, fetchTestForFunction } from '../api'
import type { TestFilePaths } from '../api'
import { useHorizontalSplit } from '../hooks/useResize'

interface Props {
  file: FileDetail | null
  fn: FunctionRecord | null
  onTestsUpdate: () => void
}

// ── helpers ───────────────────────────────────────────────────────────────────

function escapeHtml(text: string): string {
  return text
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
}

function highlightMmdb(code: string, symbols: string[]): string {
  if (!symbols.length) return escapeHtml(code)
  let result = escapeHtml(code)
  const sorted = [...symbols].sort((a, b) => b.length - a.length)
  for (const sym of sorted) {
    const esc = escapeHtml(sym)
    result = result.split(esc).join(`<mark class="mmdb-symbol">${esc}</mark>`)
  }
  return result
}

function addLineNumbers(code: string, startLine: number): string {
  return code
    .split('\n')
    .map((line, i) => {
      const n = String(startLine + i).padStart(5, ' ')
      return `<span class="select-none text-gray-600 mr-3">${n}</span>${line}`
    })
    .join('\n')
}

/**
 * Parse streamed "both" response into mmdb and gemmi sections.
 * Looks for === MMDB TEST === and === GEMMI TEST === delimiters.
 */
function parseBothSections(raw: string): { mmdb: string; gemmi: string } {
  const mmdbIdx = raw.indexOf('=== MMDB TEST ===')
  const gemmiIdx = raw.indexOf('=== GEMMI TEST ===')

  if (mmdbIdx === -1 && gemmiIdx === -1) {
    // No delimiters yet — return everything as mmdb section (still streaming)
    return { mmdb: raw, gemmi: '' }
  }
  if (mmdbIdx !== -1 && gemmiIdx === -1) {
    return { mmdb: raw.slice(mmdbIdx + '=== MMDB TEST ==='.length).trimStart(), gemmi: '' }
  }
  if (mmdbIdx === -1 && gemmiIdx !== -1) {
    return { mmdb: '', gemmi: raw.slice(gemmiIdx + '=== GEMMI TEST ==='.length).trimStart() }
  }
  const mmdbRaw = raw.slice(mmdbIdx + '=== MMDB TEST ==='.length, gemmiIdx).trim()
  const gemmiRaw = raw.slice(gemmiIdx + '=== GEMMI TEST ==='.length).trim()
  return { mmdb: mmdbRaw, gemmi: gemmiRaw }
}

function stripOuterFences(code: string): string {
  // Remove wrapping ```cpp ... ``` fences if present
  const fenced = /^```(?:cpp)?\s*\n([\s\S]*?)```\s*$/.exec(code.trim())
  return fenced ? fenced[1] : code
}

function statusColor(s: string) {
  switch (s) {
    case 'done': return 'bg-green-700 text-green-100'
    case 'reviewed': return 'bg-blue-700 text-blue-100'
    default: return 'bg-gray-600 text-gray-300'
  }
}

// ── Sub-components ────────────────────────────────────────────────────────────

function TestEditor({
  label,
  value,
  onChange,
  streaming,
  streamLabel,
  badge,
}: {
  label: string
  value: string
  onChange: (v: string) => void
  streaming?: boolean
  streamLabel?: string
  badge?: React.ReactNode
}) {
  return (
    <div className="flex flex-col overflow-hidden min-w-0 min-h-0 h-full">
      <div className="px-3 py-1.5 bg-gray-800/60 border-b border-gray-700 flex items-center gap-2 flex-shrink-0">
        <span className="text-xs font-semibold text-gray-300 uppercase tracking-wider">{label}</span>
        {badge}
        {streaming && streamLabel && (
          <span className="ml-auto text-xs text-blue-400 animate-pulse">{streamLabel}</span>
        )}
      </div>
      <div className="flex-1 relative overflow-hidden">
        <textarea
          value={value}
          onChange={e => onChange(e.target.value)}
          spellCheck={false}
          className="absolute inset-0 w-full h-full resize-none bg-gray-950 text-xs font-mono text-green-300 p-3 focus:outline-none leading-5 border-none"
          placeholder={`// ${label} code will appear here\n// You can edit it directly`}
        />
        {streaming && (
          <div className="absolute bottom-2 right-2 pointer-events-none">
            <span className="animate-pulse text-blue-400 text-lg">▋</span>
          </div>
        )}
      </div>
    </div>
  )
}

// ── Main component ────────────────────────────────────────────────────────────

export default function TestPanel({ file, fn, onTestsUpdate }: Props) {
  const [sourceCode, setSourceCode] = useState('')
  const [sourceLoading, setSourceLoading] = useState(false)
  const [showSource, setShowSource] = useState(true)

  // Resize state
  const [sourceHeight, setSourceHeight] = useState(200)
  const [terminalHeight, setTerminalHeight] = useState(180)
  const [leftPercent, editorsRef, startEditorSplitDrag] = useHorizontalSplit(50)
  // Track whether last mousedown on source handle was a drag (to avoid toggling on drag)
  const sourceDraggedRef = useRef(false)

  const [mmdbTest, setMmdbTest] = useState('')
  const [gemmiTest, setGemmiTest] = useState('')
  const [notes, setNotes] = useState('')
  const [testStatus, setTestStatus] = useState<'draft' | 'reviewed' | 'done'>('draft')

  const [model, setModel] = useState('llama3.2')
  const [extraInstructions, setExtraInstructions] = useState('')

  const [streaming, setStreaming] = useState<false | 'mmdb' | 'gemmi' | 'both'>(false)
  const [error, setError] = useState<string | null>(null)
  const [saved, setSaved] = useState(false)
  const [dirty, setDirty] = useState(false)

  // File paths & terminal
  const [filePaths, setFilePaths] = useState<TestFilePaths | null>(null)
  const [terminal, setTerminal] = useState('')
  const [compiling, setCompiling] = useState<false | 'mmdb' | 'gemmi'>(false)
  const [showTerminal, setShowTerminal] = useState(false)
  const terminalRef = useRef<HTMLDivElement>(null)

  const abortRef = useRef<AbortController | null>(null)
  const compileAbortRef = useRef<AbortController | null>(null)

  const testKey = file && fn ? `${file.rel_path}::${fn.name}:${fn.line}` : null

  // Load source when function changes
  useEffect(() => {
    if (!file || !fn) {
      setSourceCode('')
      return
    }
    setSourceLoading(true)
    fetchSource(file.rel_path, fn.line, fn.end_line)
      .then(d => setSourceCode(d.lines))
      .catch(e => setError(String(e)))
      .finally(() => setSourceLoading(false))
  }, [file, fn])

  // Load saved test when function changes — reads tests.json first, falls back to .cc files on disk
  const loadSaved = useCallback(async (relPath: string, fnName: string, fnLine: number) => {
    try {
      const rec = await fetchTestForFunction(relPath, fnName, fnLine)
      setMmdbTest(rec.mmdb_test)
      setGemmiTest(rec.gemmi_test)
      setNotes(rec.notes)
      setTestStatus(rec.status as 'draft' | 'reviewed' | 'done')
      setSaved(!!(rec.mmdb_test || rec.gemmi_test))
      setDirty(false)
    } catch (e) {
      console.error('Failed to load test', e)
    }
  }, [])

  useEffect(() => {
    if (file && fn) {
      loadSaved(file.rel_path, fn.name, fn.line)
    }
  }, [file, fn, loadSaved])

  // Load file paths when function changes
  useEffect(() => {
    if (!file || !fn) { setFilePaths(null); return }
    fetchTestFilePaths(file.rel_path, fn.name)
      .then(setFilePaths)
      .catch(() => setFilePaths(null))
  }, [file, fn])

  // Auto-scroll terminal to bottom
  useEffect(() => {
    if (terminalRef.current) {
      terminalRef.current.scrollTop = terminalRef.current.scrollHeight
    }
  }, [terminal])

  // Mark dirty on user edits
  const handleMmdbChange = (v: string) => { setMmdbTest(v); setDirty(true); setSaved(false) }
  const handleGemmiChange = (v: string) => { setGemmiTest(v); setDirty(true); setSaved(false) }
  const handleNotesChange = (v: string) => { setNotes(v); setDirty(true); setSaved(false) }

  const generate = async (target: 'mmdb' | 'gemmi' | 'both') => {
    if (!file || !fn || !sourceCode) return
    if (streaming) {
      abortRef.current?.abort()
      return
    }

    setStreaming(target)
    setError(null)
    if (target === 'mmdb' || target === 'both') setMmdbTest('')
    if (target === 'gemmi' || target === 'both') setGemmiTest('')

    abortRef.current = new AbortController()

    try {
      const response = await fetch('/api/generate-test', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          function_name: fn.name,
          source_code: sourceCode,
          mmdb_symbols: fn.mmdb_symbols,
          target,
          model,
          additional_instructions: extraInstructions,
        }),
        signal: abortRef.current.signal,
      })

      if (!response.ok) throw new Error(`Server error: ${response.status}`)
      if (!response.body) throw new Error('No response body')

      const reader = response.body.getReader()
      const decoder = new TextDecoder()
      let accumulated = ''

      while (true) {
        const { done, value } = await reader.read()
        if (done) break
        const chunk = decoder.decode(value, { stream: true })
        accumulated += chunk

        if (target === 'both') {
          const { mmdb, gemmi } = parseBothSections(accumulated)
          setMmdbTest(stripOuterFences(mmdb))
          setGemmiTest(stripOuterFences(gemmi))
        } else if (target === 'mmdb') {
          setMmdbTest(stripOuterFences(accumulated))
        } else {
          setGemmiTest(stripOuterFences(accumulated))
        }
      }

      setDirty(true)
      setSaved(false)

      // ── Auto-write to .cc files after generation completes ────────────────
      const finalMmdb = target === 'gemmi' ? '' :
        stripOuterFences(target === 'both' ? parseBothSections(accumulated).mmdb : accumulated)
      const finalGemmi = target === 'mmdb' ? '' :
        stripOuterFences(target === 'both' ? parseBothSections(accumulated).gemmi : accumulated)

      if ((finalMmdb || finalGemmi) && file && fn) {
        setShowTerminal(true)
        setTerminal(prev => prev + `\n[GENERATE COMPLETE] Writing to disk…\n`)
        try {
          const result = await writeTestFiles(
            file.rel_path, fn.name, fn.line, finalMmdb, finalGemmi,
          )
          for (const [variant, path] of Object.entries(result.written)) {
            const chars = (variant === 'mmdb' ? finalMmdb : finalGemmi).length
            setTerminal(prev => prev + `[WRITTEN]  ${path}  (${chars} chars)\n`)
          }
          // Refresh existence indicators
          const fp = await fetchTestFilePaths(file.rel_path, fn.name)
          setFilePaths(fp)
        } catch {
          setTerminal(prev => prev + `[WRITE ERROR] Could not write files — use the Write button manually\n`)
        }
      }
    } catch (e: unknown) {
      if (e instanceof Error && e.name === 'AbortError') {
        // user cancelled
      } else {
        setError(String(e))
      }
    } finally {
      setStreaming(false)
    }
  }

  const handleSave = async () => {
    if (!testKey) return
    try {
      await saveTest(testKey, mmdbTest, gemmiTest, notes, testStatus)
      setSaved(true)
      setDirty(false)
      onTestsUpdate()
    } catch (e) {
      setError(String(e))
    }
  }

  const handleStatusChange = async (s: 'draft' | 'reviewed' | 'done') => {
    setTestStatus(s)
    setDirty(true)
    setSaved(false)
  }

  // Source handle: drag to resize, click to toggle collapse
  const handleSourceHandleMouseDown = (e: React.MouseEvent) => {
    sourceDraggedRef.current = false
    const startY = e.clientY
    const startH = sourceHeight
    e.preventDefault()
    const onMove = (ev: MouseEvent) => {
      const delta = ev.clientY - startY
      if (Math.abs(delta) > 3) sourceDraggedRef.current = true
      setSourceHeight(Math.max(60, Math.min(600, startH + delta)))
      setShowSource(true)
    }
    const onUp = () => window.removeEventListener('mousemove', onMove)
    window.addEventListener('mousemove', onMove)
    window.addEventListener('mouseup', onUp, { once: true })
  }
  const handleSourceHandleClick = () => {
    if (!sourceDraggedRef.current) setShowSource(s => !s)
  }

  // Terminal handle: drag up to grow, down to shrink
  const handleTerminalHandleMouseDown = (e: React.MouseEvent) => {
    const startY = e.clientY
    const startH = terminalHeight
    e.preventDefault()
    const onMove = (ev: MouseEvent) => {
      const delta = startY - ev.clientY  // inverted
      setTerminalHeight(Math.max(60, Math.min(500, startH + delta)))
    }
    const onUp = () => window.removeEventListener('mousemove', onMove)
    window.addEventListener('mousemove', onMove)
    window.addEventListener('mouseup', onUp, { once: true })
  }

  const handleWriteFiles = async () => {
    if (!file || !fn) return
    try {
      const result = await writeTestFiles(file.rel_path, fn.name, fn.line, mmdbTest, gemmiTest)
      // Refresh file path existence flags
      const paths = await fetchTestFilePaths(file.rel_path, fn.name)
      setFilePaths(paths)
      const written = Object.values(result.written).join(', ')
      setTerminal(prev => prev + `[WRITE OK] ${written}\n`)
      setShowTerminal(true)
    } catch (e) {
      setError(String(e))
    }
  }

  const compileAndRun = async (variant: 'mmdb' | 'gemmi') => {
    if (!file || !fn) return
    const testCode = variant === 'mmdb' ? mmdbTest : gemmiTest
    if (!testCode.trim()) return

    if (compiling === variant) {
      compileAbortRef.current?.abort()
      return
    }

    setCompiling(variant)
    setTerminal('')
    setShowTerminal(true)
    setError(null)
    compileAbortRef.current = new AbortController()

    try {
      const response = await fetch('/api/tests/compile-run', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          rel_source_path: file.rel_path,
          fn_name: fn.name,
          fn_line: fn.line,
          variant,
          test_code: testCode,
        }),
        signal: compileAbortRef.current.signal,
      })

      if (!response.ok) throw new Error(`Server error: ${response.status}`)
      if (!response.body) throw new Error('No response body')

      const reader = response.body.getReader()
      const decoder = new TextDecoder()
      while (true) {
        const { done, value } = await reader.read()
        if (done) break
        setTerminal(prev => prev + decoder.decode(value, { stream: true }))
      }

      // Refresh file path flags after write-then-compile
      const paths = await fetchTestFilePaths(file.rel_path, fn.name)
      setFilePaths(paths)
    } catch (e: unknown) {
      if (!(e instanceof Error && e.name === 'AbortError')) {
        setError(String(e))
      }
    } finally {
      setCompiling(false)
    }
  }

  // ── Empty states ─────────────────────────────────────────────────────────────

  if (!file) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <div className="text-4xl mb-4">🧪</div>
        <h2 className="text-xl font-semibold text-gray-300 mb-2">Test Generator</h2>
        <p className="text-gray-500 text-sm max-w-md">
          Select a file from the sidebar, then pick a function to generate dual MMDB / Gemmi tests using a local LLM.
        </p>
      </div>
    )
  }

  if (!fn) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <div className="text-3xl mb-3">📄</div>
        <h2 className="text-lg font-semibold text-gray-300 mb-1">{file.rel_path.split('/').pop()}</h2>
        <p className="text-gray-500 text-sm">Select a function to generate tests</p>
      </div>
    )
  }

  const highlightedSource = sourceCode
    ? addLineNumbers(highlightMmdb(sourceCode, fn.mmdb_symbols), fn.line)
    : ''

  const isStreamingMmdb = streaming === 'mmdb' || streaming === 'both'
  const isStreamingGemmi = streaming === 'gemmi' || streaming === 'both'

  return (
    <div className="flex flex-col h-full overflow-hidden">
      {/* Function header */}
      <div className="px-4 py-3 bg-gray-800 border-b border-gray-700 flex-shrink-0">
        <div className="flex items-start gap-3">
          <div className="flex-1 min-w-0">
            <h3 className="text-sm font-mono font-semibold text-gray-100 truncate" title={fn.name}>
              {fn.name}
            </h3>
            <div className="flex items-center gap-3 mt-0.5 text-xs text-gray-500">
              <span>{file.rel_path}</span>
              <span>L{fn.line}–{fn.end_line}</span>
              <span>{fn.mmdb_ref_count} MMDB refs</span>
            </div>
          </div>
          <div className="flex items-center gap-2">
            {saved && !dirty && (
              <span className="text-xs text-green-400 font-medium">Saved</span>
            )}
            {dirty && (
              <span className="text-xs text-yellow-400 font-medium">Unsaved</span>
            )}
            <span className={`text-xs px-2 py-0.5 rounded font-medium ${statusColor(testStatus)}`}>
              {testStatus}
            </span>
          </div>
        </div>
        {fn.mmdb_symbols.length > 0 && (
          <div className="flex flex-wrap gap-1 mt-2">
            {fn.mmdb_symbols.map(sym => (
              <span key={sym} className="text-xs px-1.5 py-0.5 rounded bg-yellow-900/50 text-yellow-400 border border-yellow-800/50 font-mono">
                {sym}
              </span>
            ))}
          </div>
        )}
      </div>

      {/* Source code (collapsible + resizable) */}
      <div
        className="flex-shrink-0 border-b border-gray-700 overflow-auto bg-gray-950"
        style={{ height: showSource ? sourceHeight : 0, transition: showSource ? undefined : 'height 0.15s ease' }}
      >
        {sourceLoading ? (
          <div className="p-3 text-xs text-gray-500">Loading source…</div>
        ) : sourceCode ? (
          <pre
            className="p-3 text-xs font-mono leading-5 text-green-400 whitespace-pre"
            dangerouslySetInnerHTML={{ __html: highlightedSource }}
          />
        ) : (
          <div className="p-3 text-xs text-gray-500">Source not available</div>
        )}
      </div>

      {/* Source handle strip — drag to resize, click to toggle */}
      <div
        onMouseDown={handleSourceHandleMouseDown}
        onClick={handleSourceHandleClick}
        className="flex-shrink-0 flex items-center justify-center gap-1 py-0.5 border-b border-gray-700 bg-gray-800/50 hover:bg-gray-700/60 cursor-ns-resize select-none transition-colors group"
        title="Drag to resize · Click to toggle"
      >
        <span className="text-gray-600 group-hover:text-gray-400 text-xs pointer-events-none">
          {showSource ? '▲ source' : '▼ source'}
        </span>
        <span className="text-gray-700 group-hover:text-gray-500 text-xs pointer-events-none">⠿</span>
      </div>

      {/* Test editors — horizontal split with drag handle */}
      <div ref={editorsRef} className="flex flex-1 overflow-hidden min-h-0">
        <div style={{ width: `${leftPercent}%` }} className="overflow-hidden flex-shrink-0">
          <TestEditor
            label="MMDB Test"
            value={mmdbTest}
            onChange={handleMmdbChange}
            streaming={isStreamingMmdb}
            streamLabel="Generating MMDB test…"
            badge={<span className="text-xs px-1.5 py-0.5 rounded bg-yellow-900/40 text-yellow-400 border border-yellow-800/40">mmdb2</span>}
          />
        </div>

        {/* Vertical drag handle between editors */}
        <div
          onMouseDown={startEditorSplitDrag}
          className="w-1.5 flex-shrink-0 bg-gray-700 hover:bg-blue-500 active:bg-blue-400 cursor-ew-resize transition-colors select-none z-10"
          title="Drag to resize"
        />

        <div style={{ width: `${100 - leftPercent}%` }} className="overflow-hidden flex-shrink-0 border-l border-gray-700">
          <TestEditor
            label="Gemmi Test"
            value={gemmiTest}
            onChange={handleGemmiChange}
            streaming={isStreamingGemmi}
            streamLabel="Generating Gemmi test…"
            badge={<span className="text-xs px-1.5 py-0.5 rounded bg-blue-900/40 text-blue-400 border border-blue-800/40">gemmi</span>}
          />
        </div>
      </div>

      {/* File paths bar */}
      {filePaths && (
        <div className="flex-shrink-0 border-t border-gray-700 bg-gray-800/30 px-4 py-2 flex items-center gap-4 text-xs">
          <span className="text-gray-500 flex-shrink-0">Output files:</span>
          <div className="flex flex-col gap-0.5 min-w-0 flex-1">
            <div className="flex items-center gap-1.5">
              <span className={`w-1.5 h-1.5 rounded-full flex-shrink-0 ${filePaths.mmdb_exists ? 'bg-green-500' : 'bg-gray-600'}`} />
              <span className="font-mono text-gray-400 truncate" title={filePaths.mmdb}>{filePaths.mmdb}</span>
            </div>
            <div className="flex items-center gap-1.5">
              <span className={`w-1.5 h-1.5 rounded-full flex-shrink-0 ${filePaths.gemmi_exists ? 'bg-green-500' : 'bg-gray-600'}`} />
              <span className="font-mono text-gray-400 truncate" title={filePaths.gemmi}>{filePaths.gemmi}</span>
            </div>
          </div>
        </div>
      )}

      {/* Terminal output (collapsible + resizable) */}
      {showTerminal && (
        <div className="flex-shrink-0 border-t border-gray-700 flex flex-col" style={{ height: terminalHeight }}>
          {/* Terminal drag handle — drag up to grow */}
          <div
            onMouseDown={handleTerminalHandleMouseDown}
            className="h-1.5 flex-shrink-0 bg-gray-700 hover:bg-blue-500 active:bg-blue-400 cursor-ns-resize transition-colors select-none"
            title="Drag to resize terminal"
          />
          <div className="flex items-center justify-between px-3 py-1 bg-gray-900 border-b border-gray-700 flex-shrink-0">
            <span className="text-xs font-semibold text-gray-400 uppercase tracking-wider">Terminal</span>
            <div className="flex gap-2 items-center">
              {compiling && (
                <span className="text-xs text-yellow-400 animate-pulse">
                  Compiling {compiling}…
                </span>
              )}
              <button
                onClick={() => setTerminal('')}
                className="text-xs text-gray-500 hover:text-gray-300 px-1"
                title="Clear"
              >Clear</button>
              <button
                onClick={() => setShowTerminal(false)}
                className="text-xs text-gray-500 hover:text-gray-300 px-1"
                title="Hide terminal"
              >✕</button>
            </div>
          </div>
          <div
            ref={terminalRef}
            className="flex-1 overflow-auto text-xs font-mono leading-5 p-3 bg-black text-green-400 whitespace-pre-wrap"
          >
            {terminal || <span className="text-gray-600">No output yet</span>}
            {compiling && <span className="animate-pulse text-yellow-400">▋</span>}
          </div>
        </div>
      )}

      {/* Notes bar */}
      <div className="flex-shrink-0 border-t border-gray-700 bg-gray-800/40 px-4 py-2">
        <input
          type="text"
          value={notes}
          onChange={e => handleNotesChange(e.target.value)}
          placeholder="Notes (optional — describe assumptions, edge cases, etc.)"
          className="w-full bg-gray-700/60 border border-gray-600 rounded px-2 py-1.5 text-xs text-gray-200 placeholder-gray-500 focus:outline-none focus:border-blue-500"
        />
      </div>

      {/* Controls */}
      <div className="flex-shrink-0 border-t border-gray-700 bg-gray-800 p-4">
        {error && (
          <div className="mb-3 px-3 py-2 bg-red-900/30 border border-red-700/50 rounded text-xs text-red-400">
            {error}
          </div>
        )}

        <div className="flex gap-3 mb-3">
          <div className="flex flex-col gap-1">
            <label className="text-xs text-gray-400">Ollama model</label>
            <input
              type="text"
              value={model}
              onChange={e => setModel(e.target.value)}
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1.5 text-sm text-gray-100 focus:outline-none focus:border-blue-500 w-36"
            />
          </div>
          <div className="flex flex-col gap-1 flex-1">
            <label className="text-xs text-gray-400">Additional instructions</label>
            <textarea
              value={extraInstructions}
              onChange={e => setExtraInstructions(e.target.value)}
              placeholder="e.g. use TEST_F with a fixture, include edge case for empty structure…"
              rows={2}
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1.5 text-sm text-gray-100 placeholder-gray-500 focus:outline-none focus:border-blue-500 resize-none"
            />
          </div>
        </div>

        {/* Generate row */}
        <div className="flex gap-2 flex-wrap mb-2">
          <button
            onClick={() => generate('both')}
            disabled={!sourceCode}
            className="px-4 py-2 bg-blue-600 hover:bg-blue-500 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors"
          >
            {streaming === 'both' ? 'Stop' : 'Generate Both'}
          </button>
          <button
            onClick={() => generate('mmdb')}
            disabled={!sourceCode || !!streaming}
            className="px-3 py-2 bg-yellow-700 hover:bg-yellow-600 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors text-yellow-100"
          >
            {streaming === 'mmdb' ? 'Stop' : 'MMDB Only'}
          </button>
          <button
            onClick={() => generate('gemmi')}
            disabled={!sourceCode || !!streaming}
            className="px-3 py-2 bg-indigo-700 hover:bg-indigo-600 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors text-indigo-100"
          >
            {streaming === 'gemmi' ? 'Stop' : 'Gemmi Only'}
          </button>

          <div className="flex-1" />

          {/* Status selector */}
          <div className="flex items-center gap-1">
            {(['draft', 'reviewed', 'done'] as const).map(s => (
              <button
                key={s}
                onClick={() => handleStatusChange(s)}
                className={`px-2.5 py-2 rounded text-xs font-medium transition-colors border ${
                  testStatus === s
                    ? `${statusColor(s)} border-transparent`
                    : 'bg-gray-700 text-gray-400 border-gray-600 hover:bg-gray-600'
                }`}
              >
                {s}
              </button>
            ))}
          </div>

          <button
            onClick={handleSave}
            disabled={!dirty && saved}
            className="px-4 py-2 bg-green-700 hover:bg-green-600 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors"
          >
            {saved && !dirty ? '✓ Saved' : 'Save'}
          </button>
        </div>

        {/* File / compile row */}
        <div className="flex gap-2 flex-wrap">
          <button
            onClick={handleWriteFiles}
            disabled={!mmdbTest && !gemmiTest}
            className="px-3 py-1.5 bg-gray-600 hover:bg-gray-500 disabled:bg-gray-700 disabled:text-gray-500 rounded text-xs font-medium transition-colors"
          >
            Write to File
          </button>
          <button
            onClick={() => compileAndRun('mmdb')}
            disabled={!mmdbTest || compiling === 'gemmi'}
            className="px-3 py-1.5 bg-yellow-800 hover:bg-yellow-700 disabled:bg-gray-700 disabled:text-gray-500 rounded text-xs font-medium transition-colors text-yellow-100"
          >
            {compiling === 'mmdb' ? 'Stop' : 'Compile & Run MMDB'}
          </button>
          <button
            onClick={() => compileAndRun('gemmi')}
            disabled={!gemmiTest || compiling === 'mmdb'}
            className="px-3 py-1.5 bg-indigo-800 hover:bg-indigo-700 disabled:bg-gray-700 disabled:text-gray-500 rounded text-xs font-medium transition-colors text-indigo-100"
          >
            {compiling === 'gemmi' ? 'Stop' : 'Compile & Run Gemmi'}
          </button>
          {!showTerminal && terminal && (
            <button
              onClick={() => setShowTerminal(true)}
              className="px-3 py-1.5 bg-gray-700 hover:bg-gray-600 rounded text-xs font-medium transition-colors text-gray-300"
            >
              Show Output
            </button>
          )}
          <div className="flex-1" />
          <button
            onClick={() => navigator.clipboard.writeText(mmdbTest).catch(console.error)}
            disabled={!mmdbTest}
            className="px-3 py-1.5 bg-gray-600 hover:bg-gray-500 disabled:bg-gray-700 disabled:text-gray-500 rounded text-xs font-medium transition-colors"
          >
            Copy MMDB
          </button>
          <button
            onClick={() => navigator.clipboard.writeText(gemmiTest).catch(console.error)}
            disabled={!gemmiTest}
            className="px-3 py-1.5 bg-gray-600 hover:bg-gray-500 disabled:bg-gray-700 disabled:text-gray-500 rounded text-xs font-medium transition-colors"
          >
            Copy Gemmi
          </button>
        </div>
      </div>
    </div>
  )
}
