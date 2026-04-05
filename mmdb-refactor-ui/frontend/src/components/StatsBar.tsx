import type { Stats } from '../types'

interface Props {
  stats: Stats | null
}

function fmt(n: number): string {
  return n.toLocaleString()
}

export default function StatsBar({ stats }: Props) {
  return (
    <div className="flex items-center gap-6 px-6 py-3 bg-gray-800 border-b border-gray-700 text-sm flex-shrink-0 min-w-0">
      <span className="font-semibold text-blue-400 text-base whitespace-nowrap">MMDB Refactor Manager</span>
      <div className="flex-1" />
      {stats ? (
        <>
          <StatChip label="Files scanned" value={fmt(stats.total_files_scanned)} color="text-gray-300" />
          <StatChip
            label="With MMDB refs"
            value={`${fmt(stats.files_with_mmdb_refs)} (${((stats.files_with_mmdb_refs / stats.total_files_scanned) * 100).toFixed(1)}%)`}
            color="text-orange-400"
          />
          <StatChip label="Total refs" value={fmt(stats.total_mmdb_refs)} color="text-red-400" />
        </>
      ) : (
        <span className="text-gray-500 text-xs">Loading stats…</span>
      )}
      <div className="flex items-center gap-3 border-l border-gray-700 pl-4 text-xs text-gray-400">
        <span className="flex items-center gap-1"><span className="inline-block w-2 h-2 rounded-full bg-green-500" />done</span>
        <span className="flex items-center gap-1"><span className="inline-block w-2 h-2 rounded-full bg-blue-500" />in_progress</span>
        <span className="flex items-center gap-1"><span className="inline-block w-2 h-2 rounded-full bg-gray-500" />skipped</span>
      </div>
    </div>
  )
}

function StatChip({ label, value, color }: { label: string; value: string; color: string }) {
  return (
    <div className="flex flex-col items-end">
      <span className={`font-mono font-bold ${color}`}>{value}</span>
      <span className="text-gray-500 text-xs">{label}</span>
    </div>
  )
}
