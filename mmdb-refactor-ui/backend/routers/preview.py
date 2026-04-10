"""Endpoints: /api/preview-prompt/* — return rendered prompts without calling the LLM."""
from fastapi import APIRouter
from pydantic import BaseModel

from prompts import STRATEGY_CONTEXT, build_refactor, build_generate_test

router = APIRouter()


class PreviewRefactorRequest(BaseModel):
    function_name: str
    code: str
    mmdb_symbols: list[str]
    additional_instructions: str = ""


class PreviewTestRequest(BaseModel):
    function_name: str
    source_code: str
    mmdb_symbols: list[str]
    target: str = "both"   # "mmdb" | "gemmi" | "both"
    additional_instructions: str = ""
    rel_source_path: str = ""


@router.post("/api/preview-prompt/refactor")
def preview_refactor_prompt(req: PreviewRefactorRequest):
    return {
        "system": STRATEGY_CONTEXT,
        "prompt": build_refactor(req.function_name, req.code, req.mmdb_symbols,
                                 req.additional_instructions),
    }


@router.post("/api/preview-prompt/generate-test")
def preview_test_prompt(req: PreviewTestRequest):
    return {
        "system": STRATEGY_CONTEXT,
        "prompt": build_generate_test(req.function_name, req.source_code, req.mmdb_symbols,
                                      req.target, req.additional_instructions,
                                      rel_source_path=req.rel_source_path),
    }
