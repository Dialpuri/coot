import type { FileDetail, FunctionRecord, ProgressMap } from '../types'

interface Props {
  file: FileDetail | null
  onSelectFunction: (fn: FunctionRecord) => void
  selectedFn: FunctionRecord | null
  progress: ProgressMap
}

function progressBadge(status: string | undefined) {
  switch (status) {
    case 'done':
      return <span className="px-1.5 py-0.5 rounded text-xs bg-green-700 text-green-100">done</span>
    case 'in_progress':
      return <span className="px-1.5 py-0.5 rounded text-xs bg-blue-700 text-blue-100">in progress</span>
    case 'skipped':
      return <span className="px-1.5 py-0.5 rounded text-xs bg-gray-600 text-gray-300">skipped</span>
    default:
      return <span className="px-1.5 py-0.5 rounded text-xs bg-gray-700 text-gray-400">todo</span>
  }
}

function shortFnName(name: string): string {
  // Show just the last two segments of a qualified name like coot::util::foo
  const parts = name.split('::')
  if (parts.length > 2) {
    return parts.slice(-2).join('::')
  }
  return name
}

export default function FunctionList({ file, onSelectFunction, selectedFn, progress }: Props) {
  if (!file) {
    return (
      <div className="flex flex-col h-full">
        <div className="px-3 py-3 border-b border-gray-700 flex-shrink-0">
          <h2 className="text-xs font-semibold text-gray-400 uppercase tracking-wider">Functions</h2>
        </div>
        <div className="flex-1 flex items-center justify-center text-gray-600 text-sm px-4 text-center">
          Select a file from the sidebar to browse its functions
        </div>
      </div>
    )
  }

  const functions = file.functions ?? []

  return (
    <div className="flex flex-col h-full">
      <div className="px-3 py-3 border-b border-gray-700 flex-shrink-0">
        <h2 className="text-xs font-semibold text-gray-400 uppercase tracking-wider mb-1">
          Functions ({functions.length})
        </h2>
        <p className="text-xs text-gray-500 truncate font-mono" title={file.rel_path}>
          {file.rel_path.split('/').pop()}
        </p>
        <div className="mt-1 flex gap-3 text-xs text-gray-600">
          <span>{file.total_mmdb_refs} refs</span>
          <span>{file.mmdb_includes.length} includes</span>
        </div>
        {file.mmdb_includes.length > 0 && (
          <div className="mt-1.5 flex flex-wrap gap-1">
            {file.mmdb_includes.slice(0, 2).map(inc => (
              <span key={inc} className="text-xs px-1 py-0.5 rounded bg-gray-700 text-gray-400 font-mono truncate max-w-full" title={inc}>
                {inc.replace(/.*\//, '').replace(/[">]/g, '')}
              </span>
            ))}
            {file.mmdb_includes.length > 2 && (
              <span className="text-xs text-gray-600">+{file.mmdb_includes.length - 2}</span>
            )}
          </div>
        )}
      </div>

      <div className="flex-1 overflow-y-auto">
        {functions.length === 0 && (
          <div className="px-3 py-8 text-center text-xs text-gray-500">
            No functions with MMDB refs
          </div>
        )}
        {functions.map(fn => {
          const key = `${file.rel_path}::${fn.name}`
          const status = progress[key]
          const isSelected = selectedFn?.name === fn.name && selectedFn?.line === fn.line
          const displaySymbols = fn.mmdb_symbols.slice(0, 2)
          const extraSymbols = fn.mmdb_symbols.length - 2

          return (
            <button
              key={`${fn.name}:${fn.line}`}
              onClick={() => onSelectFunction(fn)}
              title={fn.name}
              className={`w-full text-left px-3 py-2.5 border-b border-gray-700/50 hover:bg-gray-700/50 transition-colors ${
                isSelected ? 'bg-blue-900/40 border-l-2 border-l-blue-500' : ''
              }`}
            >
              <div className="flex items-start gap-2">
                <div className="flex-1 min-w-0">
                  <span className="text-xs font-mono text-gray-200 block truncate">
                    {shortFnName(fn.name)}
                  </span>
                  <span className="text-xs text-gray-500">
                    L{fn.line}–{fn.end_line}
                  </span>
                </div>
                <span className="text-xs font-bold px-1.5 py-0.5 rounded bg-orange-700 text-orange-100 flex-shrink-0">
                  {fn.mmdb_ref_count}
                </span>
              </div>

              <div className="flex flex-wrap gap-1 mt-1.5">
                {displaySymbols.map(sym => (
                  <span
                    key={sym}
                    className="text-xs px-1 py-0.5 rounded bg-gray-700 text-gray-400 font-mono truncate max-w-full"
                    title={sym}
                  >
                    {sym.replace('mmdb::', '')}
                  </span>
                ))}
                {extraSymbols > 0 && (
                  <span className="text-xs px-1 py-0.5 rounded bg-gray-700 text-gray-500">
                    +{extraSymbols} more
                  </span>
                )}
              </div>

              <div className="mt-1.5">
                {progressBadge(status)}
              </div>
            </button>
          )
        })}
      </div>
    </div>
  )
}
