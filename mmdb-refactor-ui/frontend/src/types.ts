export interface FileSummary {
  rel_path: string;
  total_mmdb_refs: number;
  function_count: number;
  includes_mmdb: boolean;
}

export interface FunctionRecord {
  name: string;
  line: number;
  end_line: number;
  mmdb_ref_count: number;
  mmdb_symbols: string[];
  refs: [number, string][];
}

export interface FileDetail {
  rel_path: string;
  includes_mmdb: boolean;
  mmdb_includes: string[];
  functions: FunctionRecord[];
  file_level_refs: [number, string][];
  total_mmdb_refs: number;
}

export interface Stats {
  total_files_scanned: number;
  files_with_mmdb_refs: number;
  total_mmdb_refs: number;
}

export interface FilesResponse {
  items: FileSummary[];
  total: number;
  page: number;
  page_size: number;
}

export type ProgressStatus = 'done' | 'in_progress' | 'skipped' | 'todo';
export type ProgressMap = Record<string, ProgressStatus>;

export interface TestRecord {
  mmdb_test: string;
  gemmi_test: string;
  notes: string;
  status: 'draft' | 'reviewed' | 'done';
}
export type TestsMap = Record<string, TestRecord>;

export type AppMode = 'refactor' | 'tests';
