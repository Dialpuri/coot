import { useState, useEffect, useCallback } from 'react'
import type { FileSummary, FileDetail } from '../types'
import { fetchFiles, fetchFile } from '../api'

interface Props {
  onSelectFile: (file: FileDetail) => void
  selectedPath: string | null
}

function shortName(rel_path: string): string {
  const parts = rel_path.split('/')
  return parts.length >= 2 ? parts.slice(-2).join('/') : rel_path
}

export default function Sidebar({ onSelectFile, selectedPath }: Props) {
  const [search, setSearch] = useState('')
  const [items, setItems] = useState<FileSummary[]>([])
  const [total, setTotal] = useState(0)
  const [page, setPage] = useState(1)
  const [loading, setLoading] = useState(false)
  const [loadingFile, setLoadingFile] = useState<string | null>(null)

  const PAGE_SIZE = 50

  const load = useCallback(async (s: string, p: number, reset: boolean) => {
    setLoading(true)
    try {
      const data = await fetchFiles(s, p, PAGE_SIZE)
      setTotal(data.total)
      setItems(prev => reset ? data.items : [...prev, ...data.items])
    } catch (e) {
      console.error('Failed to load files', e)
    } finally {
      setLoading(false)
    }
  }, [])

  useEffect(() => {
    const timer = setTimeout(() => {
      setPage(1)
      load(search, 1, true)
    }, 300)
    return () => clearTimeout(timer)
  }, [search, load])

  const handleLoadMore = () => {
    const nextPage = page + 1
    setPage(nextPage)
    load(search, nextPage, false)
  }

  const handleClick = async (item: FileSummary) => {
    if (loadingFile === item.rel_path) return
    setLoadingFile(item.rel_path)
    try {
      onSelectFile(await fetchFile(item.rel_path))
    } catch (e) {
      console.error('Failed to load file detail', e)
    } finally {
      setLoadingFile(null)
    }
  }

  return (
    <div className="flex flex-col h-full bg-zinc-900">
      <div className="px-3 pt-3 pb-2 flex-shrink-0">
        <div className="flex items-baseline justify-between mb-2">
          <span className="text-xs font-medium text-zinc-400">Files</span>
          <span className="text-xs text-zinc-400">{total.toLocaleString()}</span>
        </div>
        <input
          type="text"
          placeholder="Search…"
          value={search}
          onChange={e => setSearch(e.target.value)}
          className="w-full bg-zinc-800 border border-zinc-700 rounded px-2 py-1.5 text-xs text-zinc-100 placeholder-zinc-600 focus:outline-none focus:border-blue-500"
        />
      </div>

      <div className="flex-1 overflow-y-auto">
        {items.map(item => {
          const isSelected = item.rel_path === selectedPath
          const isLoading = loadingFile === item.rel_path
          return (
            <button
              key={item.rel_path}
              onClick={() => handleClick(item)}
              title={item.rel_path}
              className={`w-full text-left px-3 py-2 border-b border-zinc-800/60 transition-colors ${
                isSelected
                  ? 'bg-zinc-800 border-l-2 border-l-blue-500 pl-[10px]'
                  : 'hover:bg-zinc-800/50'
              }`}
            >
              <div className="flex items-center gap-2">
                <span className="flex-1 text-xs text-zinc-200 truncate font-mono">
                  {isLoading ? <span className="text-blue-400">Loading…</span> : shortName(item.rel_path)}
                </span>
                <span className="text-xs text-zinc-400 tabular-nums flex-shrink-0">{item.total_mmdb_refs}</span>
              </div>
              <div className="flex items-center gap-2 mt-0.5">
                <span className="text-xs text-zinc-600 truncate flex-1">
                  {item.rel_path.split('/').slice(0, -2).join('/')}
                </span>
                <span className="text-xs text-zinc-500">{item.function_count} fn</span>
              </div>
            </button>
          )
        })}

        {loading && (
          <div className="px-3 py-4 text-center text-xs text-zinc-600">Loading…</div>
        )}

        {!loading && items.length < total && (
          <button
            onClick={handleLoadMore}
            className="w-full px-3 py-3 text-xs text-blue-500 hover:text-blue-400 transition-colors"
          >
            Load more ({total - items.length} remaining)
          </button>
        )}

        {!loading && items.length === 0 && (
          <div className="px-3 py-8 text-center text-xs text-zinc-600">No files found</div>
        )}
      </div>
    </div>
  )
}
