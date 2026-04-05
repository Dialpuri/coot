import { useState, useEffect, useCallback } from 'react'
import type { FileSummary, FileDetail } from '../types'
import { fetchFiles, fetchFile } from '../api'

interface Props {
  onSelectFile: (file: FileDetail) => void
  selectedPath: string | null
}

function refBadgeColor(count: number): string {
  if (count > 50) return 'bg-red-900 text-red-300'
  if (count > 10) return 'bg-orange-900 text-orange-300'
  return 'bg-yellow-900 text-yellow-300'
}

/** Last two path segments, e.g. "utils/foo.cc" */
function shortName(rel_path: string): string {
  const parts = rel_path.split('/')
  if (parts.length >= 2) {
    return parts.slice(-2).join('/')
  }
  return rel_path
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
      if (reset) {
        setItems(data.items)
      } else {
        setItems(prev => [...prev, ...data.items])
      }
    } catch (e) {
      console.error('Failed to load files', e)
    } finally {
      setLoading(false)
    }
  }, [])

  // Debounced search
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
      const detail = await fetchFile(item.rel_path)
      onSelectFile(detail)
    } catch (e) {
      console.error('Failed to load file detail', e)
    } finally {
      setLoadingFile(null)
    }
  }

  const hasMore = items.length < total

  return (
    <div className="flex flex-col h-full">
      {/* Header */}
      <div className="px-3 py-3 border-b border-gray-700 flex-shrink-0">
        <h2 className="text-xs font-semibold text-gray-400 uppercase tracking-wider mb-2">
          Files ({total.toLocaleString()})
        </h2>
        <input
          type="text"
          placeholder="Search files…"
          value={search}
          onChange={e => setSearch(e.target.value)}
          className="w-full bg-gray-700 border border-gray-600 rounded px-2 py-1.5 text-sm text-gray-100 placeholder-gray-500 focus:outline-none focus:border-blue-500"
        />
      </div>

      {/* List */}
      <div className="flex-1 overflow-y-auto">
        {items.map(item => {
          const isSelected = item.rel_path === selectedPath
          const isLoading = loadingFile === item.rel_path
          return (
            <button
              key={item.rel_path}
              onClick={() => handleClick(item)}
              title={item.rel_path}
              className={`w-full text-left px-3 py-2 border-b border-gray-700/50 hover:bg-gray-700/50 transition-colors ${
                isSelected ? 'bg-blue-900/40 border-l-2 border-l-blue-500' : ''
              }`}
            >
              <div className="flex items-center gap-2">
                <span className="flex-1 text-xs text-gray-200 truncate font-mono">
                  {isLoading ? (
                    <span className="text-blue-400">Loading…</span>
                  ) : (
                    shortName(item.rel_path)
                  )}
                </span>
                <span className={`text-xs font-bold px-1.5 py-0.5 rounded ${refBadgeColor(item.total_mmdb_refs)}`}>
                  {item.total_mmdb_refs}
                </span>
              </div>
              <div className="flex items-center gap-2 mt-0.5">
                <span className="text-xs text-gray-500 truncate flex-1">{item.rel_path.split('/').slice(0, -1).join('/')}</span>
                <span className="text-xs text-gray-500">{item.function_count} fn</span>
              </div>
            </button>
          )
        })}

        {loading && (
          <div className="px-3 py-4 text-center text-xs text-gray-500">Loading…</div>
        )}

        {!loading && hasMore && (
          <button
            onClick={handleLoadMore}
            className="w-full px-3 py-3 text-xs text-blue-400 hover:text-blue-300 hover:bg-gray-700/30 transition-colors"
          >
            Load more ({total - items.length} remaining)
          </button>
        )}

        {!loading && items.length === 0 && (
          <div className="px-3 py-8 text-center text-xs text-gray-500">No files found</div>
        )}
      </div>
    </div>
  )
}
