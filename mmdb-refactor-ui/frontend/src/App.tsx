import { useState, useEffect, useCallback } from 'react'
import type { FileDetail, FunctionRecord, ProgressMap, Stats, TestsMap, AppMode } from './types'
import { fetchStats, fetchProgress, fetchTests, fetchFile } from './api'
import { useHorizontalResize } from './hooks/useResize'
import StatsBar from './components/StatsBar'
import Sidebar from './components/Sidebar'
import FunctionList from './components/FunctionList'
import RefactorPanel from './components/RefactorPanel'
import TestPanel from './components/TestPanel'

const SK = 'mmdb-refactor' // localStorage key prefix

function lsGet(key: string) {
  try { return localStorage.getItem(`${SK}:${key}`) } catch { return null }
}
function lsSet(key: string, val: string) {
  try { localStorage.setItem(`${SK}:${key}`, val) } catch { /* ignore */ }
}

export default function App() {
  const [mode, setMode] = useState<AppMode>(() => (lsGet('mode') as AppMode) ?? 'refactor')
  const [stats, setStats] = useState<Stats | null>(null)
  const [sidebarWidth, startSidebarDrag] = useHorizontalResize(288, 150, 500, `${SK}:sidebar-width`)
  const [fnListWidth, startFnListDrag] = useHorizontalResize(384, 180, 600, `${SK}:fnlist-width`)
  const [selectedFile, setSelectedFile] = useState<FileDetail | null>(null)
  const [selectedFunction, setSelectedFunction] = useState<FunctionRecord | null>(null)
  const [progress, setProgress] = useState<ProgressMap>({})
  const [tests, setTests] = useState<TestsMap>({})
  const [restoring, setRestoring] = useState(true)

  const loadStats = useCallback(async () => {
    try { setStats(await fetchStats()) } catch (e) { console.error(e) }
  }, [])

  const loadProgress = useCallback(async () => {
    try { setProgress(await fetchProgress()) } catch (e) { console.error(e) }
  }, [])

  const loadTests = useCallback(async () => {
    try { setTests(await fetchTests()) } catch (e) { console.error(e) }
  }, [])

  // ── Restore selection from localStorage on mount ──────────────────────────
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
          } catch { /* ignore malformed json */ }
        }
      })
      .catch(() => { /* file may no longer be in report */ })
      .finally(() => setRestoring(false))
  }, []) // intentionally empty — runs once on mount

  useEffect(() => {
    loadStats()
    loadProgress()
    loadTests()
  }, [loadStats, loadProgress, loadTests])

  // ── Persist selection changes ──────────────────────────────────────────────
  useEffect(() => {
    if (!restoring) lsSet('mode', mode)
  }, [mode, restoring])

  useEffect(() => {
    if (!restoring && selectedFile) lsSet('selected-file', selectedFile.rel_path)
  }, [selectedFile, restoring])

  useEffect(() => {
    if (!restoring && selectedFunction) {
      lsSet('selected-fn', JSON.stringify({ name: selectedFunction.name, line: selectedFunction.line }))
    }
  }, [selectedFunction, restoring])

  const handleSelectFile = (file: FileDetail) => {
    setSelectedFile(file)
    setSelectedFunction(null)
    lsSet('selected-fn', '') // clear stale fn when file changes
  }

  // Summary counts for the mode tab badges
  const testsDone = Object.values(tests).filter(t => t.status === 'done').length
  const testsTotal = Object.keys(tests).length
  const refactorDone = Object.values(progress).filter(s => s === 'done').length

  return (
    <div className="flex flex-col h-screen bg-gray-900 text-gray-100 overflow-hidden">
      <StatsBar stats={stats} />

      {/* Mode tab bar */}
      <div className="flex items-center gap-0 px-4 bg-gray-850 border-b border-gray-700 flex-shrink-0" style={{ background: '#161b22' }}>
        <button
          onClick={() => setMode('refactor')}
          className={`px-5 py-2.5 text-sm font-medium transition-colors border-b-2 -mb-px ${
            mode === 'refactor'
              ? 'border-blue-500 text-blue-400'
              : 'border-transparent text-gray-500 hover:text-gray-300'
          }`}
        >
          Refactor
          {refactorDone > 0 && (
            <span className="ml-2 text-xs px-1.5 py-0.5 rounded-full bg-blue-900/60 text-blue-300">
              {refactorDone}
            </span>
          )}
        </button>
        <button
          onClick={() => setMode('tests')}
          className={`px-5 py-2.5 text-sm font-medium transition-colors border-b-2 -mb-px ${
            mode === 'tests'
              ? 'border-green-500 text-green-400'
              : 'border-transparent text-gray-500 hover:text-gray-300'
          }`}
        >
          Tests
          {testsTotal > 0 && (
            <span className="ml-2 text-xs px-1.5 py-0.5 rounded-full bg-green-900/60 text-green-300">
              {testsDone}/{testsTotal}
            </span>
          )}
        </button>

        {/* Test progress bar */}
        {mode === 'tests' && testsTotal > 0 && (
          <div className="flex items-center gap-2 ml-auto text-xs text-gray-500">
            <span>{testsDone} done · {testsTotal - testsDone} remaining</span>
            <div className="w-32 h-1.5 bg-gray-700 rounded-full overflow-hidden">
              <div
                className="h-full bg-green-600 rounded-full transition-all"
                style={{ width: `${(testsDone / testsTotal) * 100}%` }}
              />
            </div>
          </div>
        )}

        {restoring && (
          <span className="ml-auto text-xs text-gray-600 animate-pulse">Restoring session…</span>
        )}
      </div>

      <div className="flex flex-1 overflow-hidden">
        {/* Left column: File browser */}
        <div className="flex-shrink-0 overflow-hidden flex flex-col" style={{ width: sidebarWidth }}>
          <Sidebar
            onSelectFile={handleSelectFile}
            selectedPath={selectedFile?.rel_path ?? null}
          />
        </div>

        {/* Sidebar resize handle */}
        <div
          onMouseDown={startSidebarDrag}
          className="w-1.5 flex-shrink-0 bg-gray-700 hover:bg-blue-500 active:bg-blue-400 cursor-ew-resize transition-colors select-none z-10"
          title="Drag to resize"
        />

        {/* Middle column: Function list */}
        <div className="flex-shrink-0 overflow-hidden flex flex-col" style={{ width: fnListWidth }}>
          <FunctionList
            file={selectedFile}
            onSelectFunction={setSelectedFunction}
            selectedFn={selectedFunction}
            progress={progress}
            tests={tests}
            mode={mode}
          />
        </div>

        {/* Function list resize handle */}
        <div
          onMouseDown={startFnListDrag}
          className="w-1.5 flex-shrink-0 bg-gray-700 hover:bg-blue-500 active:bg-blue-400 cursor-ew-resize transition-colors select-none z-10"
          title="Drag to resize"
        />

        {/* Right column: Panel (mode-dependent) */}
        <div className="flex-1 overflow-hidden flex flex-col min-w-0">
          {mode === 'refactor' ? (
            <RefactorPanel
              file={selectedFile}
              fn={selectedFunction}
              onProgressUpdate={loadProgress}
              stats={stats}
            />
          ) : (
            <TestPanel
              file={selectedFile}
              fn={selectedFunction}
              onTestsUpdate={loadTests}
            />
          )}
        </div>
      </div>
    </div>
  )
}
