import type { Stats } from '../types'

interface Props {
  stats: Stats | null
}

export default function StatsBar({ stats }: Props) {
  return (
    <div className="flex items-center gap-2 px-4 py-2 bg-zinc-900 border-b border-zinc-800 text-xs text-zinc-400 flex-shrink-0 min-w-0">
      <span className="font-semibold text-zinc-100 text-sm mr-2">Coot MMDB to Gemmi Refactor Manager</span>
      {stats ? (
        <>
          <span>{stats.total_files_scanned.toLocaleString()} files</span>
          <span className="text-zinc-700">·</span>
          <span>
            {stats.files_with_mmdb_refs.toLocaleString()} with MMDB refs
            <span className="text-zinc-600 ml-1">
              ({((stats.files_with_mmdb_refs / stats.total_files_scanned) * 100).toFixed(0)}%)
            </span>
          </span>
          <span className="text-zinc-700">·</span>
          <span>{stats.total_mmdb_refs.toLocaleString()} total refs</span>
        </>
      ) : (
        <span className="text-zinc-600">Loading…</span>
      )}
    </div>
  )
}
