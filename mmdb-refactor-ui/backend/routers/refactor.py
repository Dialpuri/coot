"""Endpoint: /api/refactor (streaming Ollama refactor)"""
import json

import httpx
from fastapi import APIRouter
from fastapi.responses import StreamingResponse

from config import OLLAMA_URL
from models import RefactorRequest
from prompts import STRATEGY_CONTEXT, build_refactor

router = APIRouter()


@router.post("/api/refactor")
async def refactor(req: RefactorRequest):
    prompt = build_refactor(req.function_name, req.code, req.mmdb_symbols,
                            req.additional_instructions)
    payload = {
        "model":  req.model,
        "system": STRATEGY_CONTEXT,
        "prompt": prompt,
        "stream": True,
    }

    async def stream_ollama():
        async with httpx.AsyncClient(timeout=300.0) as client:
            async with client.stream("POST", OLLAMA_URL, json=payload) as response:
                if response.status_code != 200:
                    body = await response.aread()
                    yield f"Error from Ollama ({response.status_code}): {body.decode()}"
                    return
                async for line in response.aiter_lines():
                    if line.strip():
                        try:
                            data = json.loads(line)
                            chunk = data.get("response", "")
                            if chunk:
                                yield chunk
                        except json.JSONDecodeError:
                            pass

    return StreamingResponse(stream_ollama(), media_type="text/plain")
