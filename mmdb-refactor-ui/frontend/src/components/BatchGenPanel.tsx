import { useEffect, useRef } from 'react'
import type { BatchState, BatchConfig } from '../hooks/useBatch'

function VariantTag({ label, status, attempts }: { label: string; status: string; attempts: number }) {
  if (status === 'skip') return <span className="text-zinc-700">{label}·</span>
  const color = status === 'pass' ? 'text-emerald-500' : 'text-red-400'
  const retryColor = attempts > 1 ? 'text-amber-500' : 'text-zinc-600'
  return (
    <span className={color}>
      {label}
      {attempts > 1 && <span className={`${retryColor} text-[10px]`}>×{attempts}</span>}
    </span>
  )
}

interface Props {
  state: BatchState
  config: BatchConfig
  onConfigChange: (c: BatchConfig) => void
  onStart: () => void
  onStop: () => void
  onReset: () => void
  onClose: () => void
}

export default function BatchGenPanel({ state, config, onConfigChange, onStart, onStop, onReset, onClose }: Props) {
  const logEndRef = useRef<HTMLDivElement>(null)
  const { phase, total, done, skipped, errors, currentFn, currentFile, log } = state

  // Auto-scroll log as entries arrive
  useEffect(() => {
    logEndRef.current?.scrollIntoView({ behavior: 'smooth' })
  }, [log.length])

  const processed = done + skipped + errors
  const progress = total > 0 ? processed / total : 0

  return (
    <div className="flex flex-col h-full bg-zinc-950 overflow-hidden">
      {/* Header */}
      <div className="flex items-center justify-between px-4 py-3 border-b border-zinc-800 flex-shrink-0">
        <div>
          <h2 className="text-sm font-semibold text-zinc-100">Batch Test Generation</h2>
          <p className="text-xs text-zinc-500 mt-0.5">
            Generate tests for every function using Ollama
          </p>
        </div>
        <button onClick={onClose} className="btn btn-ghost btn-sm text-zinc-500" title="Hide panel (job continues running)">
          ×
        </button>
      </div>

      {/* Config — only when idle */}
      {phase === 'idle' && (
        <div className="px-4 py-4 border-b border-zinc-800 flex-shrink-0 flex flex-col gap-3">
          <div className="flex gap-3">
            <div className="flex flex-col gap-1">
              <label className="text-xs text-zinc-500">Model</label>
              <input
                type="text"
                value={config.model}
                onChange={e => onConfigChange({ ...config, model: e.target.value })}
                className="bg-zinc-800 border border-zinc-700 rounded px-2 py-1.5 text-sm text-zinc-100 focus:outline-none focus:border-blue-500 w-36"
              />
            </div>
            <div className="flex flex-col gap-1">
              <label className="text-xs text-zinc-500">Target</label>
              <select
                value={config.target}
                onChange={e => onConfigChange({ ...config, target: e.target.value })}
                className="bg-zinc-800 border border-zinc-700 rounded px-2 py-1.5 text-sm text-zinc-100 focus:outline-none focus:border-blue-500"
              >
                <option value="mmdb">MMDB only</option>
                <option value="gemmi">Gemmi only</option>
                <option value="both">Both</option>
              </select>
            </div>
            <div className="flex flex-col gap-1 flex-1">
              <label className="text-xs text-zinc-500">Additional instructions (applied to every function)</label>
              <input
                type="text"
                value={config.extraInstructions}
                onChange={e => onConfigChange({ ...config, extraInstructions: e.target.value })}
                placeholder="e.g. keep tests concise, use TEST_F with a shared fixture…"
                className="bg-zinc-800 border border-zinc-700 rounded px-2 py-1.5 text-sm text-zinc-100 placeholder-zinc-600 focus:outline-none focus:border-blue-500"
              />
            </div>
          </div>
          <label className="flex items-center gap-2 cursor-pointer select-none w-fit">
            <input
              type="checkbox"
              checked={config.skipExisting}
              onChange={e => onConfigChange({ ...config, skipExisting: e.target.checked })}
              className="w-3.5 h-3.5 accent-blue-500"
            />
            <span className="text-xs text-zinc-400">Skip functions that already have both tests</span>
          </label>
        </div>
      )}

      {/* Progress bar */}
      {phase !== 'idle' && (
        <div className="px-4 py-3 border-b border-zinc-800 flex-shrink-0">
          {phase === 'running' && (
            <div className="min-h-[2.5rem] mb-1.5">
              <div className="text-xs font-mono text-zinc-300 truncate">
                {currentFn ?? 'Starting…'}
              </div>
              {state.currentAttempt ? (
                <div className="flex items-start gap-1.5 mt-0.5">
                  <span className="text-xs text-amber-500 flex-shrink-0">
                    retrying {state.currentAttempt.variant} ({state.currentAttempt.attempt}/{state.currentAttempt.max})
                  </span>
                  <span className="text-xs text-zinc-600 font-mono truncate" title={state.currentAttempt.error}>
                    {state.currentAttempt.error.split('\n')[0]}
                  </span>
                </div>
              ) : (
                <div className="text-xs text-zinc-600 truncate font-mono mt-0.5">{currentFile ?? ''}</div>
              )}
            </div>
          )}
          {phase === 'finished' && (
            <div className="text-xs text-zinc-400 mb-1.5">Finished</div>
          )}
          <div className="h-1.5 bg-zinc-800 rounded-full overflow-hidden">
            <div
              className={`h-full rounded-full transition-all duration-300 ${
                phase === 'finished' && errors === 0 ? 'bg-emerald-500' : 'bg-blue-500'
              }`}
              style={{ width: `${progress * 100}%` }}
            />
          </div>
          <div className="flex gap-4 mt-2 text-xs text-zinc-600">
            <span className="tabular-nums">{processed} / {total}</span>
            <span><span className="text-emerald-500">{done}</span> generated</span>
            <span><span className="text-zinc-400">{skipped}</span> skipped</span>
            {errors > 0 && <span><span className="text-red-400">{errors}</span> errors</span>}
          </div>
        </div>
      )}

      {/* Log */}
      {log.length > 0 ? (
        <div className="flex-1 overflow-y-auto px-4 py-2 min-h-0">
          <div className="space-y-0.5">
            {log.map((entry, i) => (
              <div key={i} className="flex items-center gap-2 text-xs py-0.5">
                {entry.type === 'done' && (
                  <>
                    <span className={`w-3 flex-shrink-0 ${
                      entry.mmdb_status === 'fail' || entry.gemmi_status === 'fail' ? 'text-amber-500' : 'text-emerald-500'
                    }`}>✓</span>
                    <span className="font-mono text-zinc-300 truncate flex-1">
                      {entry.fn.split('::').pop()}
                    </span>
                    <span className="flex gap-1.5 flex-shrink-0 font-mono text-xs">
                      <VariantTag label="M" status={entry.mmdb_status ?? 'skip'} attempts={entry.mmdb_attempts ?? 0} />
                      <VariantTag label="G" status={entry.gemmi_status ?? 'skip'} attempts={entry.gemmi_attempts ?? 0} />
                    </span>
                  </>
                )}
                {entry.type === 'skip' && (
                  <>
                    <span className="text-zinc-700 w-3 flex-shrink-0">·</span>
                    <span className="font-mono text-zinc-600 truncate flex-1">
                      {entry.fn.split('::').pop()}
                    </span>
                    <span className="text-zinc-700 flex-shrink-0">skip</span>
                  </>
                )}
                {entry.type === 'error' && (
                  <>
                    <span className="text-red-500 w-3 flex-shrink-0">✗</span>
                    <span className="font-mono text-zinc-400 truncate flex-1">
                      {entry.fn ? entry.fn.split('::').pop() : 'error'}
                    </span>
                    <span className="text-red-400/70 truncate max-w-[240px] flex-shrink-0" title={entry.message}>
                      {entry.message}
                    </span>
                  </>
                )}
              </div>
            ))}
            <div ref={logEndRef} />
          </div>
        </div>
      ) : (
        <div className="flex-1 flex items-center justify-center text-zinc-700 text-sm">
          {phase === 'idle' ? 'Configure above and press Start' : 'Working…'}
        </div>
      )}

      {/* Footer */}
      <div className="flex items-center gap-2 px-4 py-3 border-t border-zinc-800 flex-shrink-0">
        {phase === 'idle' && (
          <button onClick={onStart} className="btn btn-primary">
            Start Generation
          </button>
        )}
        {phase === 'running' && (
          <button onClick={onStop} className="btn btn-secondary text-red-400 hover:text-red-300">
            Stop
          </button>
        )}
        {phase === 'finished' && (
          <>
            <span className="text-xs text-zinc-500">
              {done} generated · {skipped} skipped{errors > 0 ? ` · ${errors} errors` : ''}
            </span>
            <button onClick={onReset} className="btn btn-secondary btn-sm ml-2">
              Run Again
            </button>
          </>
        )}
        <div className="flex-1" />
        <button onClick={onClose} className="btn btn-ghost btn-sm text-zinc-500">
          {phase === 'running' ? 'Hide' : 'Close'}
        </button>
      </div>
    </div>
  )
}
