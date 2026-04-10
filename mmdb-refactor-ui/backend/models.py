"""Pydantic request/response models."""
from pydantic import BaseModel


class RefactorRequest(BaseModel):
    function_name: str
    code: str
    mmdb_symbols: list[str]
    model: str = "codellama"
    additional_instructions: str = ""


class ProgressUpdate(BaseModel):
    key: str
    status: str


class TestRecord(BaseModel):
    key: str
    mmdb_test: str = ""
    gemmi_test: str = ""
    notes: str = ""
    status: str = "draft"  # draft | reviewed | done


class GenerateTestRequest(BaseModel):
    function_name: str
    source_code: str
    mmdb_symbols: list[str]
    target: str = "both"   # "mmdb" | "gemmi" | "both"
    model: str = "gemma4"
    additional_instructions: str = ""
    rel_source_path: str = ""   # e.g. "coot-utils/coot-coord-utils.cc" — used to resolve the real header


class WriteTestRequest(BaseModel):
    rel_source_path: str   # e.g. "coot-utils/coot-coord-utils.cc"
    fn_name: str           # e.g. "coot::util::delete_residue_references_in_header_info"
    fn_line: int
    mmdb_test: str = ""
    gemmi_test: str = ""


class CompileRunRequest(BaseModel):
    rel_source_path: str
    fn_name: str
    fn_line: int
    variant: str           # "mmdb" | "gemmi"
    test_code: str         # current textarea content (written before compile)


class GitCommitRequest(BaseModel):
    rel_source_path: str
    fn_name: str
    variant: str           # "mmdb" | "gemmi" | "both"
    commit_message: str


class GenerateAllRequest(BaseModel):
    model: str = "gemma4"
    skip_existing: bool = True
    additional_instructions: str = ""


class ValidateFixRequest(BaseModel):
    rel_source_path: str
    fn_name: str
    fn_line: int
    variant: str        # "mmdb" | "gemmi"
    test_code: str
    model: str = "gemma4"
