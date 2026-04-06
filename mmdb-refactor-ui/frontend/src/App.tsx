import { useState, useEffect, useCallback } from 'react'
import type { FileDetail, FunctionRecord, ProgressMap, Stats, TestsMap, AppMode } from './types'
import { fetchStats, fetchProgress, fetchTests, fetchFile } from './api'
import { useHorizontalResize } from './hooks/useResize'
import { useBatch } from './hooks/useBatch'
import StatsBar from './components/StatsBar'
import Sidebar from './components/Sidebar'
import FunctionList from './components/FunctionList'
import RefactorPanel from './components/RefactorPanel'
import TestPanel from './components/TestPanel'
import BatchGenPanel from './components/BatchGenPanel'

const SK = 'mmdb-refactor'

function lsGet(key: string) {
  try { return localStorage.getItem(`${SK}:${key}`) } catch { return null }
}
function lsSet(key: string, val: string) {
  try { localStorage.setItem(`${SK}:${key}`, val) } catch { /* ignore */ }
}

export default function App() {
  const [mode, setMode] = useState<AppMode>(() => (lsGet('mode') as AppMode) ?? 'refactor')
  const [stats, setStats] = useState<Stats | null>(null)
  const [sidebarWidth, startSidebarDrag] = useHorizontalResize(272, 150, 480, `${SK}:sidebar-width`)
  const [fnListWidth, startFnListDrag] = useHorizontalResize(360, 180, 560, `${SK}:fnlist-width`)
  const [selectedFile, setSelectedFile] = useState<FileDetail | null>(null)
  const [selectedFunction, setSelectedFunction] = useState<FunctionRecord | null>(null)
  const [progress, setProgress] = useState<ProgressMap>({})
  const [tests, setTests] = useState<TestsMap>({})
  const [restoring, setRestoring] = useState(true)
  const [batchOpen, setBatchOpen] = useState(false)

  const loadStats = useCallback(async () => {
    try { setStats(await fetchStats()) } catch (e) { console.error(e) }
  }, [])

  const loadProgress = useCallback(async () => {
    try { setProgress(await fetchProgress()) } catch (e) { console.error(e) }
  }, [])

  const loadTests = useCallback(async () => {
    try { setTests(await fetchTests()) } catch (e) { console.error(e) }
  }, [])

  const batch = useBatch(loadTests)

  useEffect(() => {
    const savedPath = lsGet('selected-file')
    const savedFnJson = lsGet('selected-fn')
    if (!savedPath) { setRestoring(false); return }
    fetchFile(savedPath)
      .then(file => {
        setSelectedFile(file)
        if (savedFnJson) {
          try {
            const { name, line } = JSON.parse(savedFnJson) as { name: string; line: number }
            const fn = file.functions.find(f => f.name === name && f.line === line)
            if (fn) setSelectedFunction(fn)
          } catch { /* ignore */ }
        }
      })
      .catch(() => {})
      .finally(() => setRestoring(false))
  }, [])

  useEffect(() => {
    loadStats()
    loadProgress()
    loadTests()
  }, [loadStats, loadProgress, loadTests])

  useEffect(() => { if (!restoring) lsSet('mode', mode) }, [mode, restoring])
  useEffect(() => { if (!restoring && selectedFile) lsSet('selected-file', selectedFile.rel_path) }, [selectedFile, restoring])
  useEffect(() => {
    if (!restoring && selectedFunction)
      lsSet('selected-fn', JSON.stringify({ name: selectedFunction.name, line: selectedFunction.line }))
  }, [selectedFunction, restoring])

  const handleSelectFile = (file: FileDetail) => {
    setSelectedFile(file)
    setSelectedFunction(null)
    lsSet('selected-fn', '')
  }

  const testsDone = Object.values(tests).filter(t => t.status === 'done').length
  const testsTotal = Object.keys(tests).length
  const refactorDone = Object.values(progress).filter(s => s === 'done').length

  return (
    <div className="flex flex-col h-screen bg-zinc-950 text-zinc-100 overflow-hidden">
      <StatsBar stats={stats} />

      {/* Tab bar */}
      <div className="flex items-center px-2 bg-zinc-900 border-b border-zinc-800 flex-shrink-0">
        <Tab
          active={mode === 'refactor'}
          onClick={() => setMode('refactor')}
          label="Refactor"
          count={refactorDone > 0 ? `${refactorDone} done` : undefined}
        />
        <Tab
          active={mode === 'tests'}
          onClick={() => setMode('tests')}
          label="Tests"
          count={testsTotal > 0 ? `${testsDone}/${testsTotal}` : undefined}
        />

        {mode === 'tests' && (
          <button
            onClick={() => setBatchOpen(b => !b)}
            className={`ml-4 btn btn-sm ${batchOpen ? 'btn-primary' : 'btn-secondary'}`}
          >
            {batchOpen ? 'Hide batch' : 'Run All'}
            {batch.state.phase === 'running' && !batchOpen && (
              <span className="ml-1.5 flex items-center gap-1 text-blue-400">
                <span className="inline-block w-1.5 h-1.5 rounded-full bg-blue-400 animate-pulse" />
                {batch.state.done + batch.state.skipped + batch.state.errors}/{batch.state.total}
              </span>
            )}
          </button>
        )}

        {mode === 'tests' && testsTotal > 0 && (
          <div className="flex items-center gap-2 ml-auto mr-2 text-xs text-zinc-500">
            <div className="w-24 h-1 bg-zinc-800 rounded-full overflow-hidden">
              <div
                className="h-full bg-blue-600 rounded-full transition-all"
                style={{ width: `${(testsDone / testsTotal) * 100}%` }}
              />
            </div>
            <span>{testsTotal - testsDone} remaining</span>
          </div>
        )}

        {restoring && !mode && (
          <span className="ml-auto text-xs text-zinc-600 mr-2">Restoring…</span>
        )}
      </div>

      <div className="flex flex-1 overflow-hidden">
        {/* Sidebar */}
        <div className="flex-shrink-0 overflow-hidden flex flex-col border-r border-zinc-800" style={{ width: sidebarWidth }}>
          <Sidebar onSelectFile={handleSelectFile} selectedPath={selectedFile?.rel_path ?? null} />
        </div>

        <div
          onMouseDown={startSidebarDrag}
          className="w-px flex-shrink-0 bg-zinc-800 hover:bg-blue-500 cursor-ew-resize transition-colors select-none z-10"
        />

        {/* Function list */}
        <div className="flex-shrink-0 overflow-hidden flex flex-col border-r border-zinc-800" style={{ width: fnListWidth }}>
          <FunctionList
            file={selectedFile}
            onSelectFunction={setSelectedFunction}
            selectedFn={selectedFunction}
            progress={progress}
            tests={tests}
            mode={mode}
          />
        </div>

        <div
          onMouseDown={startFnListDrag}
          className="w-px flex-shrink-0 bg-zinc-800 hover:bg-blue-500 cursor-ew-resize transition-colors select-none z-10"
        />

        {/* Main panel */}
        <div className="flex-1 overflow-hidden flex flex-col min-w-0">
          {mode === 'tests' && batchOpen ? (
            <BatchGenPanel
              state={batch.state}
              config={batch.config}
              onConfigChange={batch.setConfig}
              onStart={batch.start}
              onStop={batch.stop}
              onReset={batch.reset}
              onClose={() => setBatchOpen(false)}
            />
          ) : mode === 'refactor' ? (
            <RefactorPanel file={selectedFile} fn={selectedFunction} onProgressUpdate={loadProgress} stats={stats} />
          ) : (
            <TestPanel file={selectedFile} fn={selectedFunction} onTestsUpdate={loadTests} />
          )}
        </div>
      </div>
    </div>
  )
}

function Tab({ active, onClick, label, count }: { active: boolean; onClick: () => void; label: string; count?: string }) {
  return (
    <button
      onClick={onClick}
      className={`px-4 py-2.5 text-sm font-medium border-b-2 -mb-px transition-colors ${
        active
          ? 'border-blue-500 text-zinc-100'
          : 'border-transparent text-zinc-400 hover:text-zinc-300'
      }`}
    >
      {label}
      {count && (
        <span className={`ml-2 text-xs ${active ? 'text-zinc-400' : 'text-zinc-600'}`}>{count}</span>
      )}
    </button>
  )
}
