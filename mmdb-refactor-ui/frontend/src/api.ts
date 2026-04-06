import axios from 'axios'
import type { Stats, FilesResponse, FileDetail, ProgressMap, ProgressStatus, TestsMap } from './types'

const BASE = '/api'

export async function fetchStats(): Promise<Stats> {
  const res = await axios.get<Stats>(`${BASE}/stats`)
  return res.data
}

export async function fetchFiles(
  search: string,
  page: number,
  pageSize: number
): Promise<FilesResponse> {
  const res = await axios.get<FilesResponse>(`${BASE}/files`, {
    params: { search, page, page_size: pageSize },
  })
  return res.data
}

export async function fetchFile(relPath: string): Promise<FileDetail> {
  const res = await axios.get<FileDetail>(`${BASE}/file/${relPath}`)
  return res.data
}

export async function fetchSource(
  path: string,
  start: number,
  end: number
): Promise<{ lines: string; start: number; end: number }> {
  const res = await axios.get<{ lines: string; start: number; end: number }>(
    `${BASE}/source`,
    { params: { path, start, end } }
  )
  return res.data
}

export async function fetchProgress(): Promise<ProgressMap> {
  const res = await axios.get<ProgressMap>(`${BASE}/progress`)
  return res.data
}

export async function postProgress(key: string, status: ProgressStatus | 'todo'): Promise<void> {
  await axios.post(`${BASE}/progress`, { key, status })
}

export async function fetchTests(): Promise<TestsMap> {
  const res = await axios.get<TestsMap>(`${BASE}/tests`)
  return res.data
}

export async function fetchTestForFunction(
  rel_source_path: string,
  fn_name: string,
  fn_line: number,
): Promise<{ mmdb_test: string; gemmi_test: string; notes: string; status: string }> {
  const res = await axios.get(`${BASE}/tests/load`, {
    params: { rel_source_path, fn_name, fn_line },
  })
  return res.data
}

export async function saveTest(
  key: string,
  mmdb_test: string,
  gemmi_test: string,
  notes: string,
  status: string
): Promise<void> {
  await axios.post(`${BASE}/tests`, { key, mmdb_test, gemmi_test, notes, status })
}

export interface TestFilePaths {
  mmdb: string
  gemmi: string
  mmdb_exists: boolean
  gemmi_exists: boolean
}

export async function fetchTestFilePaths(
  rel_source_path: string,
  fn_name: string
): Promise<TestFilePaths> {
  const res = await axios.get<TestFilePaths>(`${BASE}/tests/file-paths`, {
    params: { rel_source_path, fn_name },
  })
  return res.data
}

export async function writeTestFiles(
  rel_source_path: string,
  fn_name: string,
  fn_line: number,
  mmdb_test: string,
  gemmi_test: string
): Promise<{ written: Record<string, string> }> {
  const res = await axios.post(`${BASE}/tests/write`, {
    rel_source_path, fn_name, fn_line, mmdb_test, gemmi_test,
  })
  return res.data
}

export async function gitCommitTestFiles(
  rel_source_path: string,
  fn_name: string,
  variant: 'mmdb' | 'gemmi' | 'both',
  commit_message: string,
  onChunk: (chunk: string) => void,
  signal?: AbortSignal,
): Promise<void> {
  const response = await fetch(`${BASE}/tests/git-commit`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ rel_source_path, fn_name, variant, commit_message }),
    signal,
  })
  if (!response.ok) throw new Error(`Server error: ${response.status}`)
  if (!response.body) throw new Error('No response body')
  const reader = response.body.getReader()
  const decoder = new TextDecoder()
  while (true) {
    const { done, value } = await reader.read()
    if (done) break
    onChunk(decoder.decode(value, { stream: true }))
  }
}
