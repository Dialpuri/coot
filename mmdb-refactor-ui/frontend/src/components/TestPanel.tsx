import { useState, useEffect, useRef, useCallback } from 'react'
import type { FileDetail, FunctionRecord } from '../types'
import { fetchSource, saveTest, fetchTestFilePaths, writeTestFiles, fetchTestForFunction } from '../api'
import type { TestFilePaths } from '../api'
import { useHorizontalSplit } from '../hooks/useResize'
import { highlightCppWithLines, highlightForEditor } from '../highlight'
import Editor from 'react-simple-code-editor'

interface Props {
  file: FileDetail | null
  fn: FunctionRecord | null
  onTestsUpdate: () => void
}

function parseBothSections(raw: string): { mmdb: string; gemmi: string } {
  const mmdbIdx = raw.indexOf('=== MMDB TEST ===')
  const gemmiIdx = raw.indexOf('=== GEMMI TEST ===')
  if (mmdbIdx === -1 && gemmiIdx === -1) return { mmdb: raw, gemmi: '' }
  if (mmdbIdx !== -1 && gemmiIdx === -1) return { mmdb: raw.slice(mmdbIdx + '=== MMDB TEST ==='.length).trimStart(), gemmi: '' }
  if (mmdbIdx === -1 && gemmiIdx !== -1) return { mmdb: '', gemmi: raw.slice(gemmiIdx + '=== GEMMI TEST ==='.length).trimStart() }
  return {
    mmdb: raw.slice(mmdbIdx + '=== MMDB TEST ==='.length, gemmiIdx).trim(),
    gemmi: raw.slice(gemmiIdx + '=== GEMMI TEST ==='.length).trim(),
  }
}

function stripOuterFences(code: string): string {
  const fenced = /^```(?:cpp)?\s*\n([\s\S]*?)```\s*$/.exec(code.trim())
  return fenced ? fenced[1] : code
}

function TestEditor({
  label,
  value,
  onChange,
  streaming,
  streamLabel,
  tag,
}: {
  label: string
  value: string
  onChange: (v: string) => void
  streaming?: boolean
  streamLabel?: string
  tag?: string
}) {
  return (
    <div className="flex flex-col overflow-hidden min-w-0 min-h-0 h-full">
      <div className="px-3 py-1.5 bg-zinc-900/60 border-b border-zinc-800 flex items-center gap-2 flex-shrink-0">
        <span className="text-xs text-zinc-400">{label}</span>
        {tag && <span className="text-xs font-mono text-zinc-600">{tag}</span>}
        {streaming && streamLabel && (
          <span className="ml-auto text-xs text-blue-400 animate-pulse">{streamLabel}</span>
        )}
      </div>
      <div className="flex-1 relative overflow-auto bg-[#0d1117]">
        <Editor
          value={value}
          onValueChange={onChange}
          highlight={highlightForEditor}
          padding={12}
          className="rsce-container min-h-full"
          textareaClassName="focus:outline-none"
          placeholder={`// ${label} will appear here\n// You can edit it directly`}
          style={{ minHeight: '100%' }}
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

export default function TestPanel({ file, fn, onTestsUpdate }: Props) {
  const [sourceCode, setSourceCode] = useState('')
  const [sourceLoading, setSourceLoading] = useState(false)
  const [showSource, setShowSource] = useState(true)
  const [sourceHeight, setSourceHeight] = useState(200)
  const [terminalHeight, setTerminalHeight] = useState(180)
  const [leftPercent, editorsRef, startEditorSplitDrag] = useHorizontalSplit(50)
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

  const [filePaths, setFilePaths] = useState<TestFilePaths | null>(null)
  const [terminal, setTerminal] = useState('')
  const [compiling, setCompiling] = useState<false | 'mmdb' | 'gemmi'>(false)
  const [showTerminal, setShowTerminal] = useState(false)
  const terminalRef = useRef<HTMLDivElement>(null)
  const abortRef = useRef<AbortController | null>(null)
  const compileAbortRef = useRef<AbortController | null>(null)

  const testKey = file && fn ? `${file.rel_path}::${fn.name}:${fn.line}` : null

  useEffect(() => {
    if (!file || !fn) { setSourceCode(''); return }
    setSourceLoading(true)
    fetchSource(file.rel_path, fn.line, fn.end_line)
      .then(d => setSourceCode(d.lines))
      .catch(e => setError(String(e)))
      .finally(() => setSourceLoading(false))
  }, [file, fn])

  const loadSaved = useCallback(async (relPath: string, fnName: string, fnLine: number) => {
    try {
      const rec = await fetchTestForFunction(relPath, fnName, fnLine)
      setMmdbTest(rec.mmdb_test)
      setGemmiTest(rec.gemmi_test)
      setNotes(rec.notes)
      setTestStatus(rec.status as 'draft' | 'reviewed' | 'done')
      setSaved(!!(rec.mmdb_test || rec.gemmi_test))
      setDirty(false)
    } catch (e) { console.error('Failed to load test', e) }
  }, [])

  useEffect(() => {
    if (file && fn) loadSaved(file.rel_path, fn.name, fn.line)
  }, [file, fn, loadSaved])

  useEffect(() => {
    if (!file || !fn) { setFilePaths(null); return }
    fetchTestFilePaths(file.rel_path, fn.name).then(setFilePaths).catch(() => setFilePaths(null))
  }, [file, fn])

  useEffect(() => {
    if (terminalRef.current) terminalRef.current.scrollTop = terminalRef.current.scrollHeight
  }, [terminal])

  const handleMmdbChange = (v: string) => { setMmdbTest(v); setDirty(true); setSaved(false) }
  const handleGemmiChange = (v: string) => { setGemmiTest(v); setDirty(true); setSaved(false) }
  const handleNotesChange = (v: string) => { setNotes(v); setDirty(true); setSaved(false) }

  const generate = async (target: 'mmdb' | 'gemmi' | 'both') => {
    if (!file || !fn || !sourceCode) return
    if (streaming) { abortRef.current?.abort(); return }
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
          function_name: fn.name, source_code: sourceCode, mmdb_symbols: fn.mmdb_symbols,
          target, model, additional_instructions: extraInstructions,
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
        accumulated += decoder.decode(value, { stream: true })
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

      const finalMmdb = target === 'gemmi' ? '' : stripOuterFences(target === 'both' ? parseBothSections(accumulated).mmdb : accumulated)
      const finalGemmi = target === 'mmdb' ? '' : stripOuterFences(target === 'both' ? parseBothSections(accumulated).gemmi : accumulated)

      if ((finalMmdb || finalGemmi) && file && fn) {
        setShowTerminal(true)
        setTerminal(prev => prev + '\n[GENERATE COMPLETE] Writing to disk…\n')
        try {
          const result = await writeTestFiles(file.rel_path, fn.name, fn.line, finalMmdb, finalGemmi)
          for (const [variant, path] of Object.entries(result.written)) {
            const chars = (variant === 'mmdb' ? finalMmdb : finalGemmi).length
            setTerminal(prev => prev + `[WRITTEN]  ${path}  (${chars} chars)\n`)
          }
          const fp = await fetchTestFilePaths(file.rel_path, fn.name)
          setFilePaths(fp)
        } catch {
          setTerminal(prev => prev + '[WRITE ERROR] Could not write files\n')
        }
      }
    } catch (e: unknown) {
      if (!(e instanceof Error && e.name === 'AbortError')) setError(String(e))
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
    } catch (e) { setError(String(e)) }
  }

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

  const handleTerminalHandleMouseDown = (e: React.MouseEvent) => {
    const startY = e.clientY
    const startH = terminalHeight
    e.preventDefault()
    const onMove = (ev: MouseEvent) => {
      const delta = startY - ev.clientY
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
      const paths = await fetchTestFilePaths(file.rel_path, fn.name)
      setFilePaths(paths)
      setTerminal(prev => prev + `[WRITE OK] ${Object.values(result.written).join(', ')}\n`)
      setShowTerminal(true)
    } catch (e) { setError(String(e)) }
  }

  const compileAndRun = async (variant: 'mmdb' | 'gemmi') => {
    if (!file || !fn) return
    const testCode = variant === 'mmdb' ? mmdbTest : gemmiTest
    if (!testCode.trim()) return
    if (compiling === variant) { compileAbortRef.current?.abort(); return }

    setCompiling(variant)
    setTerminal('')
    setShowTerminal(true)
    setError(null)
    compileAbortRef.current = new AbortController()

    try {
      const response = await fetch('/api/tests/compile-run', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ rel_source_path: file.rel_path, fn_name: fn.name, fn_line: fn.line, variant, test_code: testCode }),
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
      const paths = await fetchTestFilePaths(file.rel_path, fn.name)
      setFilePaths(paths)
    } catch (e: unknown) {
      if (!(e instanceof Error && e.name === 'AbortError')) setError(String(e))
    } finally {
      setCompiling(false)
    }
  }

  // ── Empty states ──────────────────────────────────────────────────────────

  if (!file) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <h2 className="text-base font-semibold text-zinc-300 mb-2">Test Generator</h2>
        <p className="text-zinc-600 text-sm max-w-sm">
          Select a file and function to generate dual MMDB / Gemmi tests using a local LLM.
        </p>
      </div>
    )
  }

  if (!fn) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <h2 className="text-sm font-mono font-semibold text-zinc-300 mb-1">{file.rel_path.split('/').pop()}</h2>
        <p className="text-zinc-600 text-sm">Select a function to generate tests</p>
      </div>
    )
  }

  const highlightedSource = sourceCode ? highlightCppWithLines(sourceCode, fn.line, fn.mmdb_symbols) : ''
  const isStreamingMmdb = streaming === 'mmdb' || streaming === 'both'
  const isStreamingGemmi = streaming === 'gemmi' || streaming === 'both'

  return (
    <div className="flex flex-col h-full overflow-hidden bg-zinc-950">
      {/* Function header */}
      <div className="px-4 py-3 bg-zinc-900 border-b border-zinc-800 flex-shrink-0">
        <div className="flex items-start gap-3">
          <div className="flex-1 min-w-0">
            <h3 className="text-sm font-mono font-medium text-zinc-100 truncate" title={fn.name}>{fn.name}</h3>
            <div className="flex items-center gap-3 mt-0.5 text-xs text-zinc-600">
              <span>{file.rel_path}</span>
              <span>L{fn.line}–{fn.end_line}</span>
              <span>{fn.mmdb_ref_count} refs</span>
            </div>
          </div>
          <div className="flex items-center gap-2 text-xs">
            {saved && !dirty && <span className="text-emerald-500">saved</span>}
            {dirty && <span className="text-amber-500">unsaved</span>}
            <span className={`${
              testStatus === 'done' ? 'text-emerald-500' : testStatus === 'reviewed' ? 'text-blue-400' : 'text-zinc-500'
            }`}>{testStatus}</span>
          </div>
        </div>
        {fn.mmdb_symbols.length > 0 && (
          <div className="flex flex-wrap gap-2 mt-2">
            {fn.mmdb_symbols.map(sym => (
              <span key={sym} className="text-xs font-mono text-amber-500/80 px-1.5 py-0.5 rounded  bg-yellow-900/50 text-yellow-400 border border-yellow-800/50 font-mono">{sym}</span>
            ))}
          </div>
        )}
      </div>

      {/* Source (collapsible + resizable) */}
      <div
        className="flex-shrink-0 border-b border-zinc-800 overflow-auto"
        style={{ height: showSource ? sourceHeight : 0, transition: showSource ? undefined : 'height 0.15s ease' }}
      >
        {sourceLoading ? (
          <div className="p-3 text-xs text-zinc-600">Loading…</div>
        ) : sourceCode ? (
          <pre
            className="language-cpp p-3 text-xs leading-5 bg-[#0d1117] min-h-full whitespace-pre"
            dangerouslySetInnerHTML={{ __html: highlightedSource }}
          />
        ) : (
          <div className="p-3 text-xs text-zinc-600">Source not available</div>
        )}
      </div>

      {/* Source drag/toggle handle */}
      <div
        onMouseDown={handleSourceHandleMouseDown}
        onClick={() => { if (!sourceDraggedRef.current) setShowSource(s => !s) }}
        className="flex-shrink-0 flex items-center justify-center gap-1 py-0.5 border-b border-zinc-800 bg-zinc-900/50 hover:bg-zinc-800/60 cursor-ns-resize select-none transition-colors group"
        title="Drag to resize · Click to toggle"
      >
        <span className="text-zinc-700 group-hover:text-zinc-500 text-xs pointer-events-none">
          {showSource ? '▲' : '▼'} source
        </span>
      </div>

      {/* Editors — horizontal split */}
      <div ref={editorsRef} className="flex flex-1 overflow-hidden min-h-0">
        <div style={{ width: `${leftPercent}%` }} className="overflow-hidden flex-shrink-0">
          <TestEditor
            label="MMDB Test"
            value={mmdbTest}
            onChange={handleMmdbChange}
            streaming={isStreamingMmdb}
            streamLabel="Generating…"
            tag="mmdb2"
          />
        </div>

        <div
          onMouseDown={startEditorSplitDrag}
          className="w-px flex-shrink-0 bg-zinc-800 hover:bg-blue-500 cursor-ew-resize transition-colors select-none z-10"
        />

        <div style={{ width: `${100 - leftPercent}%` }} className="overflow-hidden flex-shrink-0">
          <TestEditor
            label="Gemmi Test"
            value={gemmiTest}
            onChange={handleGemmiChange}
            streaming={isStreamingGemmi}
            streamLabel="Generating…"
            tag="gemmi"
          />
        </div>
      </div>

      {/* File paths */}
      {filePaths && (
        <div className="flex-shrink-0 border-t border-zinc-800 bg-zinc-900/30 px-4 py-2 flex items-center gap-3 text-xs">
          <span className="text-zinc-600 flex-shrink-0">Output</span>
          <div className="flex flex-col gap-0.5 min-w-0 flex-1">
            <div className="flex items-center gap-1.5">
              <span className={`w-1.5 h-1.5 rounded-full flex-shrink-0 ${filePaths.mmdb_exists ? 'bg-emerald-500' : 'bg-zinc-700'}`} />
              <span className="font-mono text-zinc-500 truncate" title={filePaths.mmdb}>{filePaths.mmdb}</span>
            </div>
            <div className="flex items-center gap-1.5">
              <span className={`w-1.5 h-1.5 rounded-full flex-shrink-0 ${filePaths.gemmi_exists ? 'bg-emerald-500' : 'bg-zinc-700'}`} />
              <span className="font-mono text-zinc-500 truncate" title={filePaths.gemmi}>{filePaths.gemmi}</span>
            </div>
          </div>
        </div>
      )}

      {/* Terminal */}
      {showTerminal && (
        <div className="flex-shrink-0 border-t border-zinc-800 flex flex-col" style={{ height: terminalHeight }}>
          <div
            onMouseDown={handleTerminalHandleMouseDown}
            className="h-1 flex-shrink-0 bg-zinc-800 hover:bg-blue-500 cursor-ns-resize transition-colors select-none"
          />
          <div className="flex items-center justify-between px-3 py-1 bg-zinc-900 border-b border-zinc-800 flex-shrink-0">
            <span className="text-xs text-zinc-500">
              Terminal
              {compiling && <span className="ml-2 text-amber-500 animate-pulse">compiling {compiling}…</span>}
            </span>
            <div className="flex gap-3">
              <button onClick={() => setTerminal('')} className="btn btn-ghost btn-sm">Clear</button>
              <button onClick={() => setShowTerminal(false)} className="btn btn-ghost btn-sm">×</button>
            </div>
          </div>
          <div
            ref={terminalRef}
            className="flex-1 overflow-auto text-xs font-mono leading-5 p-3 bg-black text-emerald-400 whitespace-pre-wrap"
          >
            {terminal || <span className="text-zinc-700">No output yet</span>}
            {compiling && <span className="animate-pulse text-amber-400">▋</span>}
          </div>
        </div>
      )}

      {/* Notes */}
      <div className="flex-shrink-0 border-t border-zinc-800 bg-zinc-900/40 px-4 py-2">
        <input
          type="text"
          value={notes}
          onChange={e => handleNotesChange(e.target.value)}
          placeholder="Notes — assumptions, edge cases, etc."
          className="w-full bg-zinc-800/60 border border-zinc-700 rounded px-2 py-1.5 text-xs text-zinc-200 placeholder-zinc-600 focus:outline-none focus:border-blue-500"
        />
      </div>

      {/* Controls */}
      <div className="flex-shrink-0 border-t border-zinc-800 bg-zinc-900 p-4">
        {error && (
          <div className="mb-3 px-3 py-2 bg-red-950/50 border border-red-900/50 rounded text-xs text-red-400">{error}</div>
        )}

        <div className="flex gap-3 mb-3">
          <div className="flex flex-col gap-1">
            <label className="text-xs text-zinc-500">Model</label>
            <input
              type="text"
              value={model}
              onChange={e => setModel(e.target.value)}
              className="bg-zinc-800 border border-zinc-700 rounded px-2 py-1.5 text-sm text-zinc-100 focus:outline-none focus:border-blue-500 w-36"
            />
          </div>
          <div className="flex flex-col gap-1 flex-1">
            <label className="text-xs text-zinc-500">Additional instructions</label>
            <textarea
              value={extraInstructions}
              onChange={e => setExtraInstructions(e.target.value)}
              placeholder="e.g. use TEST_F with a fixture, include edge case for empty structure…"
              rows={2}
              className="bg-zinc-800 border border-zinc-700 rounded px-2 py-1.5 text-sm text-zinc-100 placeholder-zinc-600 focus:outline-none focus:border-blue-500 resize-none"
            />
          </div>
        </div>

        {/* Generate row */}
        <div className="flex gap-2 flex-wrap mb-2">
          <button onClick={() => generate('both')} disabled={!sourceCode} className="btn btn-primary">
            {streaming === 'both' ? 'Stop' : 'Generate Both'}
          </button>
          <button onClick={() => generate('mmdb')} disabled={!sourceCode || !!streaming} className="btn btn-secondary">
            {streaming === 'mmdb' ? 'Stop' : 'MMDB Only'}
          </button>
          <button onClick={() => generate('gemmi')} disabled={!sourceCode || !!streaming} className="btn btn-secondary">
            {streaming === 'gemmi' ? 'Stop' : 'Gemmi Only'}
          </button>

          <div className="flex-1" />

          {/* Status */}
          <div className="flex items-center border border-zinc-700 rounded overflow-hidden">
            {(['draft', 'reviewed', 'done'] as const).map(s => (
              <button
                key={s}
                onClick={() => { setTestStatus(s); setDirty(true); setSaved(false) }}
                className={`px-2.5 py-1.5 text-xs transition-colors ${
                  testStatus === s
                    ? s === 'done' ? 'bg-zinc-700 text-emerald-400' : s === 'reviewed' ? 'bg-zinc-700 text-blue-400' : 'bg-zinc-700 text-zinc-200'
                    : 'text-zinc-600 hover:text-zinc-400'
                }`}
              >
                {s}
              </button>
            ))}
          </div>

          <button onClick={handleSave} disabled={!dirty && saved} className="btn btn-secondary text-emerald-400 hover:text-emerald-300">
            {saved && !dirty ? 'Saved' : 'Save'}
          </button>
        </div>

        {/* File/compile row */}
        <div className="flex gap-2 flex-wrap">
          <button onClick={handleWriteFiles} disabled={!mmdbTest && !gemmiTest} className="btn btn-secondary btn-sm">
            Write to File
          </button>
          <button onClick={() => compileAndRun('mmdb')} disabled={!mmdbTest || compiling === 'gemmi'} className="btn btn-secondary btn-sm">
            {compiling === 'mmdb' ? 'Stop' : 'Run MMDB'}
          </button>
          <button onClick={() => compileAndRun('gemmi')} disabled={!gemmiTest || compiling === 'mmdb'} className="btn btn-secondary btn-sm">
            {compiling === 'gemmi' ? 'Stop' : 'Run Gemmi'}
          </button>
          {!showTerminal && terminal && (
            <button onClick={() => setShowTerminal(true)} className="btn btn-ghost btn-sm">Show output</button>
          )}
          <div className="flex-1" />
          <button onClick={() => navigator.clipboard.writeText(mmdbTest).catch(console.error)} disabled={!mmdbTest} className="btn btn-ghost btn-sm">
            Copy MMDB
          </button>
          <button onClick={() => navigator.clipboard.writeText(gemmiTest).catch(console.error)} disabled={!gemmiTest} className="btn btn-ghost btn-sm">
            Copy Gemmi
          </button>
        </div>
      </div>
    </div>
  )
}
