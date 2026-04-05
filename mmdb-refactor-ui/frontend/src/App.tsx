import { useState, useEffect, useCallback } from 'react'
import type { FileDetail, FunctionRecord, ProgressMap, Stats } from './types'
import { fetchStats, fetchProgress } from './api'
import StatsBar from './components/StatsBar'
import Sidebar from './components/Sidebar'
import FunctionList from './components/FunctionList'
import RefactorPanel from './components/RefactorPanel'

export default function App() {
  const [stats, setStats] = useState<Stats | null>(null)
  const [selectedFile, setSelectedFile] = useState<FileDetail | null>(null)
  const [selectedFunction, setSelectedFunction] = useState<FunctionRecord | null>(null)
  const [progress, setProgress] = useState<ProgressMap>({})

  const loadStats = useCallback(async () => {
    try {
      const s = await fetchStats()
      setStats(s)
    } catch (e) {
      console.error('Failed to load stats', e)
    }
  }, [])

  const loadProgress = useCallback(async () => {
    try {
      const p = await fetchProgress()
      setProgress(p)
    } catch (e) {
      console.error('Failed to load progress', e)
    }
  }, [])

  useEffect(() => {
    loadStats()
    loadProgress()
  }, [loadStats, loadProgress])

  const handleSelectFile = (file: FileDetail) => {
    setSelectedFile(file)
    setSelectedFunction(null)
  }

  const handleSelectFunction = (fn: FunctionRecord) => {
    setSelectedFunction(fn)
  }

  const handleProgressUpdate = () => {
    loadProgress()
  }

  return (
    <div className="flex flex-col h-screen bg-gray-900 text-gray-100 overflow-hidden">
      <StatsBar stats={stats} />
      <div className="flex flex-1 overflow-hidden">
        {/* Left column: File browser */}
        <div className="w-72 flex-shrink-0 border-r border-gray-700 overflow-hidden flex flex-col">
          <Sidebar
            onSelectFile={handleSelectFile}
            selectedPath={selectedFile?.rel_path ?? null}
          />
        </div>

        {/* Middle column: Function list */}
        <div className="w-96 flex-shrink-0 border-r border-gray-700 overflow-hidden flex flex-col">
          <FunctionList
            file={selectedFile}
            onSelectFunction={handleSelectFunction}
            selectedFn={selectedFunction}
            progress={progress}
          />
        </div>

        {/* Right column: Refactor panel */}
        <div className="flex-1 overflow-hidden flex flex-col">
          <RefactorPanel
            file={selectedFile}
            fn={selectedFunction}
            onProgressUpdate={handleProgressUpdate}
            stats={stats}
          />
        </div>
      </div>
    </div>
  )
}
