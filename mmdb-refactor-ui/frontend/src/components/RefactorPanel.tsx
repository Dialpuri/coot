import { useState, useEffect, useRef } from 'react'
import type { FileDetail, FunctionRecord, Stats } from '../types'
import { fetchSource, postProgress, previewRefactorPrompt } from '../api'
import type { PromptPreview } from '../api'
import { highlightCppWithLines } from '../highlight'
import PromptModal from './PromptModal'

interface Props {
  file: FileDetail | null
  fn: FunctionRecord | null
  onProgressUpdate: () => void
  stats?: Stats | null
}

export default function RefactorPanel({ file, fn, onProgressUpdate, stats }: Props) {
  const [sourceCode, setSourceCode] = useState<string>('')
  const [gemmiCode, setGemmiCode] = useState<string>('')
  const [model, setModel] = useState('gemma4')
  const [additionalInstructions, setAdditionalInstructions] = useState('')
  const [streaming, setStreaming] = useState(false)
  const [sourceLoading, setSourceLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const abortRef = useRef<AbortController | null>(null)
  const [promptPreview, setPromptPreview] = useState<PromptPreview | null>(null)
  const [promptLoading, setPromptLoading] = useState(false)
  const [promptError, setPromptError] = useState<string | null>(null)

  const showPrompt = async () => {
    if (!fn) return
    setPromptPreview(null)
    setPromptError(null)
    setPromptLoading(true)
    try {
      const data = await previewRefactorPrompt(fn.name, sourceCode, fn.mmdb_symbols, additionalInstructions)
      setPromptPreview(data)
    } catch (e) {
      setPromptError(String(e))
    } finally {
      setPromptLoading(false)
    }
  }

  useEffect(() => {
    if (!file || !fn) { setSourceCode(''); setGemmiCode(''); return }
    setSourceLoading(true)
    setError(null)
    fetchSource(file.rel_path, fn.line, fn.end_line)
      .then(data => setSourceCode(data.lines))
      .catch(e => setError(String(e)))
      .finally(() => setSourceLoading(false))
  }, [file, fn])

  const startRefactor = async (clear: boolean) => {
    if (!file || !fn) return
    if (streaming) { abortRef.current?.abort(); return }
    if (clear) setGemmiCode('')
    setStreaming(true)
    setError(null)
    abortRef.current = new AbortController()

    try {
      const response = await fetch('/api/refactor', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          function_name: fn.name,
          code: sourceCode,
          mmdb_symbols: fn.mmdb_symbols,
          model,
          additional_instructions: additionalInstructions,
        }),
        signal: abortRef.current.signal,
      })

      if (!response.ok) throw new Error(`Server error: ${response.status}`)
      if (!response.body) throw new Error('No response body')

      const reader = response.body.getReader()
      const decoder = new TextDecoder()
      let accumulated = clear ? '' : gemmiCode

      while (true) {
        const { done, value } = await reader.read()
        if (done) break
        accumulated += decoder.decode(value, { stream: true })
        setGemmiCode(accumulated)
      }
    } catch (e: unknown) {
      if (!(e instanceof Error && e.name === 'AbortError')) setError(String(e))
    } finally {
      setStreaming(false)
    }
  }

  const markProgress = async (status: 'done' | 'skipped' | 'in_progress') => {
    if (!file || !fn) return
    await postProgress(`${file.rel_path}::${fn.name}`, status)
    onProgressUpdate()
  }

  const clearProgress = async () => {
    if (!file || !fn) return
    await postProgress(`${file.rel_path}::${fn.name}`, 'todo')
    onProgressUpdate()
  }

  // ── Empty states ──────────────────────────────────────────────────────────

  if (!file) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <h2 className="text-base font-semibold text-zinc-300 mb-2">MMDB → Gemmi Refactor</h2>
        <p className="text-zinc-600 text-sm mb-6 max-w-sm">
          Select a file and function to generate a Gemmi equivalent using a local LLM.
        </p>
        {stats && (
          <div className="flex gap-6 text-sm text-zinc-500">
            <span><span className="text-zinc-300 font-mono">{stats.total_files_scanned.toLocaleString()}</span> files</span>
            <span><span className="text-zinc-300 font-mono">{stats.files_with_mmdb_refs.toLocaleString()}</span> with refs</span>
            <span><span className="text-zinc-300 font-mono">{stats.total_mmdb_refs.toLocaleString()}</span> total refs</span>
          </div>
        )}
      </div>
    )
  }

  if (!fn) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <h2 className="text-sm font-mono font-semibold text-zinc-300 mb-1">{file.rel_path.split('/').pop()}</h2>
        <p className="text-zinc-600 text-xs mb-3">{file.rel_path}</p>
        <p className="text-zinc-600 text-sm">Select a function to begin refactoring</p>
        {file.mmdb_includes.length > 0 && (
          <div className="mt-4 text-xs text-zinc-600 font-mono">
            {file.mmdb_includes.map(inc => <div key={inc}>{inc}</div>)}
          </div>
        )}
      </div>
    )
  }

  const highlightedSource = sourceCode ? highlightCppWithLines(sourceCode, fn.line, fn.mmdb_symbols) : ''

  return (
    <div className="flex flex-col h-full overflow-hidden bg-zinc-950">
      {(promptPreview !== null || promptLoading || promptError !== null) && (
        <PromptModal
          title="Prompt Preview — Refactor"
          data={promptPreview}
          loading={promptLoading}
          error={promptError}
          onClose={() => { setPromptPreview(null); setPromptError(null) }}
        />
      )}
      {/* Function header */}
      <div className="px-4 py-3 bg-zinc-900 border-b border-zinc-800 flex-shrink-0">
        <h3 className="text-sm font-mono font-medium text-zinc-100 truncate" title={fn.name}>{fn.name}</h3>
        <div className="flex items-center gap-3 mt-0.5 text-xs text-zinc-600">
          <span>{file.rel_path}</span>
          <span>L{fn.line}–{fn.end_line}</span>
          <span>{fn.mmdb_ref_count} refs</span>
        </div>
        {fn.mmdb_symbols.length > 0 && (
          <div className="flex flex-wrap gap-2 mt-2">
            {fn.mmdb_symbols.map(sym => (
              <span key={sym} className="text-xs font-mono text-amber-500/80">{sym}</span>
            ))}
          </div>
        )}
      </div>

      {/* Code panels */}
      <div className="flex flex-1 overflow-hidden">
        {/* Before */}
        <div className="flex-1 flex flex-col border-r border-zinc-800 overflow-hidden">
          <div className="px-3 py-1.5 bg-zinc-900/60 border-b border-zinc-800 flex-shrink-0">
            <span className="text-xs text-zinc-500">Before — MMDB</span>
          </div>
          <div className="flex-1 overflow-auto">
            {sourceLoading ? (
              <div className="p-4 text-xs text-zinc-600">Loading…</div>
            ) : sourceCode ? (
              <pre
                className="language-cpp p-3 text-xs leading-5 bg-[#0d1117] min-h-full whitespace-pre"
                dangerouslySetInnerHTML={{ __html: highlightedSource }}
              />
            ) : (
              <div className="p-4 text-xs text-zinc-600">Source not available</div>
            )}
          </div>
        </div>

        {/* After */}
        <div className="flex-1 flex flex-col overflow-hidden">
          <div className="px-3 py-1.5 bg-zinc-900/60 border-b border-zinc-800 flex-shrink-0 flex items-center gap-2">
            <span className="text-xs text-zinc-500">After — Gemmi</span>
            {streaming && <span className="text-xs text-blue-400 animate-pulse">Generating…</span>}
          </div>
          <div className="flex-1 overflow-auto">
            {gemmiCode ? (
              <pre className="p-3 text-xs font-mono leading-5 text-blue-300/80 bg-zinc-950 min-h-full whitespace-pre">
                {gemmiCode}
                {streaming && <span className="animate-pulse text-blue-400">▋</span>}
              </pre>
            ) : (
              <div className="flex items-center justify-center h-full text-zinc-600 text-xs">
                Click Generate to produce a Gemmi equivalent
              </div>
            )}
          </div>
        </div>
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
              value={additionalInstructions}
              onChange={e => setAdditionalInstructions(e.target.value)}
              placeholder="e.g. use const references, preserve existing comments…"
              rows={2}
              className="bg-zinc-800 border border-zinc-700 rounded px-2 py-1.5 text-sm text-zinc-100 placeholder-zinc-600 focus:outline-none focus:border-blue-500 resize-none"
            />
          </div>
        </div>

        <div className="flex gap-2 flex-wrap">
          <button
            onClick={() => startRefactor(true)}
            disabled={!sourceCode}
            className="btn btn-primary"
          >
            {streaming ? 'Stop' : 'Generate'}
          </button>
          <button
            onClick={() => startRefactor(false)}
            disabled={!sourceCode || streaming}
            className="btn btn-secondary"
          >
            Append
          </button>
          <button
            onClick={() => navigator.clipboard.writeText(gemmiCode).catch(console.error)}
            disabled={!gemmiCode}
            className="btn btn-secondary"
          >
            Copy
          </button>
          <button
            onClick={showPrompt}
            disabled={!sourceCode}
            className="btn btn-secondary"
            title="Preview the prompt that will be sent to the LLM"
          >
            Preview Prompt
          </button>
          <div className="flex-1" />
          <button onClick={() => markProgress('done')} disabled={!file || !fn} className="btn btn-secondary text-emerald-400 hover:text-emerald-300">
            Done
          </button>
          <button onClick={() => markProgress('skipped')} disabled={!file || !fn} className="btn btn-secondary">
            Skip
          </button>
          <button onClick={clearProgress} disabled={!file || !fn} className="btn btn-ghost btn-sm text-zinc-600 hover:text-zinc-400" title="Clear status">
            ×
          </button>
        </div>
      </div>
    </div>
  )
}
