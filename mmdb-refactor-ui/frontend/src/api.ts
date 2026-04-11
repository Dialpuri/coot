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

// ── Probe files ──────────────────────────────────────────────────────────────

export interface ProbeFiles {
  source: string
  prompt: string
  source_exists: boolean
  prompt_exists: boolean
}

export async function fetchProbeFiles(): Promise<ProbeFiles> {
  const res = await axios.get<ProbeFiles>(`${BASE}/probe/files`)
  return res.data
}

// ── Prompt preview ────────────────────────────────────────────────────────────

export interface PromptPreview {
  system: string
  prompt: string
}

export async function previewRefactorPrompt(
  function_name: string,
  code: string,
  mmdb_symbols: string[],
  additional_instructions: string,
): Promise<PromptPreview> {
  const res = await axios.post<PromptPreview>(`${BASE}/preview-prompt/refactor`, {
    function_name, code, mmdb_symbols, additional_instructions,
  })
  return res.data
}

export async function previewTestPrompt(
  function_name: string,
  source_code: string,
  mmdb_symbols: string[],
  target: string,
  additional_instructions: string,
  rel_source_path = '',
): Promise<PromptPreview> {
  const res = await axios.post<PromptPreview>(`${BASE}/preview-prompt/generate-test`, {
    function_name, source_code, mmdb_symbols, target, additional_instructions, rel_source_path,
  })
  return res.data
}

// ── Batch generation ──────────────────────────────────────────────────────────

export type BatchEvent =
  | { type: 'start';    total: number }
  | { type: 'progress'; key: string; fn: string; file: string; done: number; skipped: number; errors: number; total: number }
  | { type: 'skip';     key: string; fn: string; file: string; done: number; skipped: number; errors: number; total: number }
  | { type: 'attempt';  key: string; fn: string; file: string; variant: string; attempt: number; max: number; error: string; done: number; skipped: number; errors: number; total: number }
  | { type: 'done';     key: string; fn: string; file: string; done: number; skipped: number; errors: number; total: number; mmdb_status: string; mmdb_attempts: number; gemmi_status: string; gemmi_attempts: number }
  | { type: 'error';    key: string; fn: string; file: string; done: number; skipped: number; errors: number; total: number; message: string }
  | { type: 'finish';   done: number; skipped: number; errors: number; total: number }

export async function* generateAllTests(
  model: string,
  skipExisting: boolean,
  additionalInstructions: string,
  signal?: AbortSignal,
): AsyncGenerator<BatchEvent> {
  const response = await fetch(`${BASE}/tests/generate-all`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ model, skip_existing: skipExisting, additional_instructions: additionalInstructions }),
    signal,
  })
  if (!response.ok) throw new Error(`Server error: ${response.status}`)
  if (!response.body) throw new Error('No response body')

  const reader = response.body.getReader()
  const decoder = new TextDecoder()
  let buffer = ''

  while (true) {
    const { done, value } = await reader.read()
    if (done) break
    buffer += decoder.decode(value, { stream: true })
    const lines = buffer.split('\n')
    buffer = lines.pop() ?? ''
    for (const line of lines) {
      if (line.trim()) {
        try { yield JSON.parse(line) as BatchEvent } catch { /* skip malformed */ }
      }
    }
  }
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

// ── Single-function validate + fix ────────────────────────────────────────────

export type ValidateEvent =
  | { type: 'start';          max: number; variant: string }
  | { type: 'compiling';      attempt: number; max: number }
  | { type: 'compile_output'; text: string; ok: boolean }
  | { type: 'running';        attempt: number }
  | { type: 'run_output';     text: string; ok: boolean }
  | { type: 'fixing';         attempt: number }
  | { type: 'fixed_code';     code: string }
  | { type: 'done';           status: 'pass' | 'fail'; attempts: number; code: string }

export async function* validateAndFixTest(
  rel_source_path: string,
  fn_name: string,
  fn_line: number,
  variant: 'mmdb' | 'gemmi',
  test_code: string,
  model: string,
  signal?: AbortSignal,
): AsyncGenerator<ValidateEvent> {
  const response = await fetch(`${BASE}/tests/validate-fix`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ rel_source_path, fn_name, fn_line, variant, test_code, model }),
    signal,
  })
  if (!response.ok) throw new Error(`Server error: ${response.status}`)
  if (!response.body) throw new Error('No response body')

  const reader = response.body.getReader()
  const decoder = new TextDecoder()
  let buffer = ''

  while (true) {
    const { done, value } = await reader.read()
    if (done) break
    buffer += decoder.decode(value, { stream: true })
    const lines = buffer.split('\n')
    buffer = lines.pop() ?? ''
    for (const line of lines) {
      if (line.trim()) {
        try { yield JSON.parse(line) as ValidateEvent } catch { /* skip */ }
      }
    }
  }
}
