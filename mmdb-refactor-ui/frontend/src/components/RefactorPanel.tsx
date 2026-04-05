import { useState, useEffect, useRef } from 'react'
import type { FileDetail, FunctionRecord, Stats } from '../types'
import { fetchSource, postProgress } from '../api'

interface Props {
  file: FileDetail | null
  fn: FunctionRecord | null
  onProgressUpdate: () => void
  stats?: Stats | null
}

function highlightMmdbSymbols(code: string, symbols: string[]): string {
  if (!symbols.length) return escapeHtml(code)
  let result = escapeHtml(code)
  // Sort by length descending to replace longer tokens first
  const sorted = [...symbols].sort((a, b) => b.length - a.length)
  for (const sym of sorted) {
    const escaped = escapeHtml(sym)
    // Replace whole-word occurrences
    result = result.split(escaped).join(
      `<mark class="mmdb-symbol">${escaped}</mark>`
    )
  }
  return result
}

function escapeHtml(text: string): string {
  return text
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
}

function addLineNumbers(code: string, startLine: number): string {
  return code
    .split('\n')
    .map((line, i) => {
      const lineNo = String(startLine + i).padStart(5, ' ')
      return `<span class="select-none text-gray-600 mr-3">${lineNo}</span>${line}`
    })
    .join('\n')
}

export default function RefactorPanel({ file, fn, onProgressUpdate, stats }: Props) {
  const [sourceCode, setSourceCode] = useState<string>('')
  const [gemmiCode, setGemmiCode] = useState<string>('')
  const [model, setModel] = useState('codellama')
  const [additionalInstructions, setAdditionalInstructions] = useState('')
  const [streaming, setStreaming] = useState(false)
  const [sourceLoading, setSourceLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const abortRef = useRef<AbortController | null>(null)

  // Load source code when function changes
  useEffect(() => {
    if (!file || !fn) {
      setSourceCode('')
      setGemmiCode('')
      return
    }
    setSourceLoading(true)
    setError(null)
    fetchSource(file.rel_path, fn.line, fn.end_line + 100)
      .then(data => setSourceCode(data.lines))
      .catch(e => setError(String(e)))
      .finally(() => setSourceLoading(false))
  }, [file, fn])

  const startRefactor = async (clear: boolean) => {
    if (!file || !fn) return
    if (streaming) {
      abortRef.current?.abort()
      return
    }
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

      if (!response.ok) {
        throw new Error(`Server error: ${response.status}`)
      }

      if (!response.body) {
        throw new Error('No response body')
      }

      const reader = response.body.getReader()
      const decoder = new TextDecoder()

      let accumulated = clear ? '' : gemmiCode

      while (true) {
        const { done, value } = await reader.read()
        if (done) break
        const chunk = decoder.decode(value, { stream: true })
        accumulated += chunk
        setGemmiCode(accumulated)
      }
    } catch (e: unknown) {
      if (e instanceof Error && e.name === 'AbortError') {
        // User cancelled
      } else {
        setError(String(e))
      }
    } finally {
      setStreaming(false)
    }
  }

  const copyGemmiCode = () => {
    navigator.clipboard.writeText(gemmiCode).catch(console.error)
  }

  const markProgress = async (status: 'done' | 'skipped' | 'in_progress') => {
    if (!file || !fn) return
    const key = `${file.rel_path}::${fn.name}`
    await postProgress(key, status)
    onProgressUpdate()
  }

  const clearProgress = async () => {
    if (!file || !fn) return
    const key = `${file.rel_path}::${fn.name}`
    await postProgress(key, 'todo')
    onProgressUpdate()
  }

  // ── Empty states ───────────────────────────────────────────────────────────

  if (!file) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <div className="text-4xl mb-4">⚡</div>
        <h2 className="text-xl font-semibold text-gray-300 mb-2">MMDB → Gemmi Refactor Manager</h2>
        <p className="text-gray-500 text-sm mb-6 max-w-md">
          Browse files in the sidebar, select a file to see its functions, then pick a function to generate a Gemmi refactor using a local LLM.
        </p>
        {stats && (
          <div className="flex gap-6 text-sm">
            <div className="text-center">
              <div className="text-2xl font-bold text-blue-400">{stats.total_files_scanned.toLocaleString()}</div>
              <div className="text-gray-500">files scanned</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-orange-400">{stats.files_with_mmdb_refs.toLocaleString()}</div>
              <div className="text-gray-500">with MMDB refs</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-red-400">{stats.total_mmdb_refs.toLocaleString()}</div>
              <div className="text-gray-500">total refs</div>
            </div>
          </div>
        )}
      </div>
    )
  }

  if (!fn) {
    return (
      <div className="flex flex-col h-full items-center justify-center px-8 text-center">
        <div className="text-3xl mb-3">📄</div>
        <h2 className="text-lg font-semibold text-gray-300 mb-1">{file.rel_path.split('/').pop()}</h2>
        <p className="text-gray-500 text-sm mb-4">{file.rel_path}</p>
        <p className="text-gray-600 text-sm">Select a function from the middle panel to begin refactoring</p>
        <div className="mt-4 text-xs text-gray-600">
          {file.functions.length} function{file.functions.length !== 1 ? 's' : ''} · {file.total_mmdb_refs} MMDB refs
        </div>
        {file.mmdb_includes.length > 0 && (
          <div className="mt-4 max-w-md text-left">
            <div className="text-xs text-gray-500 mb-1">MMDB includes:</div>
            {file.mmdb_includes.map(inc => (
              <div key={inc} className="text-xs font-mono text-gray-400">{inc}</div>
            ))}
          </div>
        )}
      </div>
    )
  }

  // ── Main view ──────────────────────────────────────────────────────────────

  const highlightedSource = sourceCode
    ? addLineNumbers(
        highlightMmdbSymbols(sourceCode, fn.mmdb_symbols),
        fn.line
      )
    : ''

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

      {/* Code panels */}
      <div className="flex flex-1 overflow-hidden">
        {/* Before panel */}
        <div className="flex-1 flex flex-col border-r border-gray-700 overflow-hidden">
          <div className="px-3 py-1.5 bg-gray-800/50 border-b border-gray-700 flex-shrink-0">
            <span className="text-xs font-semibold text-gray-400 uppercase tracking-wider">Before (MMDB)</span>
          </div>
          <div className="flex-1 overflow-auto">
            {sourceLoading ? (
              <div className="p-4 text-xs text-gray-500">Loading source…</div>
            ) : sourceCode ? (
              <pre
                className="p-3 text-xs font-mono leading-5 text-green-400 bg-gray-950 min-h-full whitespace-pre"
                dangerouslySetInnerHTML={{ __html: highlightedSource }}
              />
            ) : (
              <div className="p-4 text-xs text-gray-500">Source not available</div>
            )}
          </div>
        </div>

        {/* After panel */}
        <div className="flex-1 flex flex-col overflow-hidden">
          <div className="px-3 py-1.5 bg-gray-800/50 border-b border-gray-700 flex-shrink-0">
            <span className="text-xs font-semibold text-gray-400 uppercase tracking-wider">After (Gemmi)</span>
            {streaming && (
              <span className="ml-2 text-xs text-blue-400 animate-pulse">Generating…</span>
            )}
          </div>
          <div className="flex-1 overflow-auto">
            {gemmiCode ? (
              <pre className="p-3 text-xs font-mono leading-5 text-blue-300 bg-gray-950 min-h-full whitespace-pre">
                {gemmiCode}
                {streaming && <span className="animate-pulse text-blue-400">▋</span>}
              </pre>
            ) : (
              <div className="flex flex-col items-center justify-center h-full text-gray-600 text-sm">
                <div className="text-2xl mb-2">→</div>
                <p>Click &lsquo;Generate&rsquo; to get Gemmi equivalent</p>
              </div>
            )}
          </div>
        </div>
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
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1.5 text-sm text-gray-100 focus:outline-none focus:border-blue-500 w-40"
            />
          </div>
          <div className="flex flex-col gap-1 flex-1">
            <label className="text-xs text-gray-400">Additional instructions (optional)</label>
            <textarea
              value={additionalInstructions}
              onChange={e => setAdditionalInstructions(e.target.value)}
              placeholder="e.g. use const references, preserve existing comments…"
              rows={2}
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1.5 text-sm text-gray-100 placeholder-gray-500 focus:outline-none focus:border-blue-500 resize-none"
            />
          </div>
        </div>

        <div className="flex gap-2 flex-wrap">
          <button
            onClick={() => startRefactor(true)}
            disabled={!sourceCode}
            className="px-4 py-2 bg-blue-600 hover:bg-blue-500 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors"
          >
            {streaming ? 'Stop' : 'Generate'}
          </button>
          <button
            onClick={() => startRefactor(false)}
            disabled={!sourceCode || streaming}
            className="px-4 py-2 bg-gray-600 hover:bg-gray-500 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors"
          >
            Regenerate
          </button>
          <button
            onClick={copyGemmiCode}
            disabled={!gemmiCode}
            className="px-4 py-2 bg-gray-600 hover:bg-gray-500 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors"
          >
            Copy
          </button>
          <div className="flex-1" />
          <button
            onClick={() => markProgress('done')}
            disabled={!file || !fn}
            className="px-4 py-2 bg-green-700 hover:bg-green-600 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors"
          >
            ✓ Done
          </button>
          <button
            onClick={() => markProgress('skipped')}
            disabled={!file || !fn}
            className="px-4 py-2 bg-yellow-700 hover:bg-yellow-600 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors"
          >
            → Skip
          </button>
          <button
            onClick={clearProgress}
            disabled={!file || !fn}
            className="px-3 py-2 bg-red-900 hover:bg-red-800 disabled:bg-gray-700 disabled:text-gray-500 rounded text-sm font-medium transition-colors text-red-300"
            title="Clear status"
          >
            ×
          </button>
        </div>
      </div>
    </div>
  )
}
