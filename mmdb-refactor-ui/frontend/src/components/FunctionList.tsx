import React from 'react'
import type { FileDetail, FunctionRecord, ProgressMap, TestsMap, AppMode } from '../types'

interface Props {
  file: FileDetail | null
  onSelectFunction: (fn: FunctionRecord) => void
  selectedFn: FunctionRecord | null
  progress: ProgressMap
  tests: TestsMap
  mode: AppMode
}

function StatusDot({ color }: { color: string }) {
  return <span className={`inline-block w-1.5 h-1.5 rounded-full flex-shrink-0 ${color}`} />
}

function progressIndicator(status: string | undefined): React.ReactNode {
  switch (status) {
    case 'done':
      return <span className="flex items-center gap-1.5 text-xs text-emerald-500"><StatusDot color="bg-emerald-500" />done</span>
    case 'in_progress':
      return <span className="flex items-center gap-1.5 text-xs text-blue-400"><StatusDot color="bg-blue-400" />in progress</span>
    case 'skipped':
      return <span className="flex items-center gap-1.5 text-xs text-zinc-500"><StatusDot color="bg-zinc-500" />skipped</span>
    default:
      return null
  }
}

function testIndicator(status: string | undefined, hasMmdb: boolean, hasGemmi: boolean): React.ReactNode {
  if (!hasMmdb && !hasGemmi) return null
  const parts: string[] = []
  if (hasMmdb) parts.push('M')
  if (hasGemmi) parts.push('G')
  const dotColor = status === 'done' ? 'bg-emerald-500' : status === 'reviewed' ? 'bg-blue-400' : 'bg-zinc-500'
  const textColor = status === 'done' ? 'text-emerald-500' : status === 'reviewed' ? 'text-blue-400' : 'text-zinc-400'
  return (
    <span className={`flex items-center gap-1.5 text-xs ${textColor}`}>
      <StatusDot color={dotColor} />
      {parts.join(' | ')} · {status ?? 'draft'}
    </span>
  )
}

function shortFnName(name: string): string {
  const parts = name.split('::')
  return parts.length > 2 ? parts.slice(-2).join('::') : name
}

export default function FunctionList({ file, onSelectFunction, selectedFn, progress, tests, mode }: Props) {
  if (!file) {
    return (
      <div className="flex flex-col h-full bg-zinc-900">
        <div className="px-3 pt-3 pb-2 border-b border-zinc-800 flex-shrink-0">
          <span className="text-xs font-medium text-zinc-400">Functions</span>
        </div>
        <div className="flex-1 flex items-center justify-center text-zinc-600 text-xs px-4 text-center">
          Select a file to browse its functions
        </div>
      </div>
    )
  }

  const functions = file.functions ?? []

  return (
    <div className="flex flex-col h-full bg-zinc-900">
      <div className="px-3 pt-3 pb-2 border-b border-zinc-800 flex-shrink-0">
        <div className="flex items-baseline justify-between">
          <span className="text-xs font-medium text-zinc-400">Functions</span>
          <span className="text-xs text-zinc-400">{functions.length}</span>
        </div>
        <p className="text-xs text-zinc-500 truncate font-mono mt-1" title={file.rel_path}>
          {file.rel_path.split('/').pop()}
        </p>
      </div>

      <div className="flex-1 overflow-y-auto">
        {functions.length === 0 && (
          <div className="px-3 py-8 text-center text-xs text-zinc-600">No functions with MMDB refs</div>
        )}
        {functions.map(fn => {
          const refactorKey = `${file.rel_path}::${fn.name}`
          const testKey = `${file.rel_path}::${fn.name}:${fn.line}`
          const refactorStatus = progress[refactorKey]
          const testRecord = tests[testKey]
          const isSelected = selectedFn?.name === fn.name && selectedFn?.line === fn.line

          return (
            <button
              key={`${fn.name}:${fn.line}`}
              onClick={() => onSelectFunction(fn)}
              title={fn.name}
              className={`w-full text-left px-3 py-2.5 border-b border-zinc-800/60 transition-colors ${
                isSelected
                  ? 'bg-zinc-800 border-l-2 border-l-blue-500 pl-[10px]'
                  : 'hover:bg-zinc-800/50'
              }`}
            >
              <div className="flex items-start gap-2">
                <div className="flex-1 min-w-0">
                  <span className="text-xs font-mono text-zinc-200 block truncate">
                    {shortFnName(fn.name)}
                  </span>
                  <span className="text-xs text-zinc-400">
                    L{fn.line}–{fn.end_line} · {fn.mmdb_ref_count} refs
                  </span>
                </div>
              </div>

              {fn.mmdb_symbols.length > 0 && (
                <div className="flex flex-wrap gap-1 mt-1.5">
                  {fn.mmdb_symbols.slice(0, 3).map(sym => (
                    <span key={sym} className="text-xs font-mono text-zinc-500 truncate max-w-full" title={sym}>
                      {sym.replace('mmdb::', '')}
                    </span>
                  ))}
                  {fn.mmdb_symbols.length > 3 && (
                    <span className="text-xs text-zinc-700">+{fn.mmdb_symbols.length - 3}</span>
                  )}
                </div>
              )}

              <div className="mt-1.5">
                {mode === 'tests'
                  ? testIndicator(testRecord?.status, !!(testRecord?.mmdb_test), !!(testRecord?.gemmi_test))
                  : progressIndicator(refactorStatus)
                }
              </div>
            </button>
          )
        })}
      </div>
    </div>
  )
}
