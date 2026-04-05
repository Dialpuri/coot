import axios from 'axios'
import type { Stats, FilesResponse, FileDetail, ProgressMap, ProgressStatus } from './types'

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
