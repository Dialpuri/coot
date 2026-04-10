"""
Ollama HTTP client.

Two entry points:

* `call_ollama`   — wait for the full response, return the assistant's text.
* `stream_ollama` — yield chunks as they arrive, tagged "thinking" or "response"
                    so callers can show the model's reasoning separately from
                    the final answer in the UI.

Both honour the global `LLM_THINK` flag (CLI: --no-think) and pass `think=true`
to Ollama when it is enabled. Models that don't support thinking simply ignore
the field.
"""
import json
from typing import AsyncIterator

import httpx

from config import LLM_THINK, OLLAMA_URL
from prompts import STRATEGY_CONTEXT


def _build_payload(model: str, prompt: str, system: str | None, *, stream: bool) -> dict:
    payload: dict = {
        "model":  model,
        "system": STRATEGY_CONTEXT if system is None else system,
        "prompt": prompt,
        "stream": stream,
    }
    if LLM_THINK:
        payload["think"] = True
    return payload


async def call_ollama(model: str, prompt: str, system: str | None = None) -> str:
    """Non-streaming call. Returns just the assistant's `response` text.

    Any `thinking` content from a reasoning model is dropped on the floor — use
    `stream_ollama` if you want to surface it.
    """
    payload = _build_payload(model, prompt, system, stream=False)
    async with httpx.AsyncClient(timeout=600.0) as client:
        resp = await client.post(OLLAMA_URL, json=payload)
        resp.raise_for_status()
        return resp.json().get("response", "")


async def stream_ollama(
    model: str,
    prompt: str,
    system: str | None = None,
) -> AsyncIterator[dict]:
    """Streaming call. Yields dicts of the form:

        {"kind": "thinking" | "response", "text": str}

    Reasoning models emit "thinking" chunks first, then "response" chunks.
    Non-thinking models only emit "response" chunks.
    """
    payload = _build_payload(model, prompt, system, stream=True)
    async with httpx.AsyncClient(timeout=600.0) as client:
        async with client.stream("POST", OLLAMA_URL, json=payload) as resp:
            if resp.status_code != 200:
                body = await resp.aread()
                raise RuntimeError(
                    f"Ollama {resp.status_code}: {body.decode(errors='replace')[:500]}"
                )
            async for line in resp.aiter_lines():
                if not line.strip():
                    continue
                try:
                    data = json.loads(line)
                except json.JSONDecodeError:
                    continue
                thinking = data.get("thinking", "")
                if thinking:
                    yield {"kind": "thinking", "text": thinking}
                response = data.get("response", "")
                if response:
                    yield {"kind": "response", "text": response}
