"""Ollama client helpers (non-streaming and streaming)."""
import json
from typing import AsyncIterator

import httpx

from config import OLLAMA_URL
from prompts import STRATEGY_CONTEXT


async def call_ollama(model: str, prompt: str, system: str | None = None) -> str:
    """Non-streaming Ollama call — waits for full response.

    `system` defaults to the MMDB→Gemmi refactor STRATEGY_CONTEXT. Callers that
    need a different framing (e.g. the probe, which must stay in MMDB-only
    mode) should pass an explicit system prompt.
    """
    payload = {
        "model":  model,
        "system": STRATEGY_CONTEXT if system is None else system,
        "prompt": prompt,
        "stream": False,
    }
    async with httpx.AsyncClient(timeout=600.0) as client:
        resp = await client.post(OLLAMA_URL, json=payload)
        resp.raise_for_status()
        return resp.json().get("response", "")


async def stream_ollama(
    model: str,
    prompt: str,
    system: str | None = None,
) -> AsyncIterator[str]:
    """Streaming Ollama call — yields text chunks as they arrive.

    Same semantics as `call_ollama`, but returns an async generator so the
    caller can show progress while the LLM is still generating. Compute on
    the Ollama side is identical to the non-streaming variant; only the
    transport differs.
    """
    payload = {
        "model":  model,
        "system": STRATEGY_CONTEXT if system is None else system,
        "prompt": prompt,
        "stream": True,
    }
    async with httpx.AsyncClient(timeout=600.0) as client:
        async with client.stream("POST", OLLAMA_URL, json=payload) as resp:
            if resp.status_code != 200:
                body = await resp.aread()
                raise RuntimeError(f"Ollama {resp.status_code}: {body.decode(errors='replace')[:500]}")
            async for line in resp.aiter_lines():
                if not line.strip():
                    continue
                try:
                    data = json.loads(line)
                except json.JSONDecodeError:
                    continue
                chunk = data.get("response", "")
                if chunk:
                    yield chunk
