"""Non-streaming Ollama client helper."""
import httpx

from config import OLLAMA_URL
from prompts import STRATEGY_CONTEXT


async def call_ollama(model: str, prompt: str) -> str:
    """Non-streaming Ollama call — waits for full response."""
    payload = {"model": model, "system": STRATEGY_CONTEXT, "prompt": prompt, "stream": False}
    async with httpx.AsyncClient(timeout=600.0) as client:
        resp = await client.post(OLLAMA_URL, json=payload)
        resp.raise_for_status()
        return resp.json().get("response", "")
