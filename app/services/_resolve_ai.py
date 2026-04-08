"""
AI 제공자 추상화 레이어

Gemini / OpenAI 를 동일한 인터페이스로 호출합니다.
acmg.py, 미래의 report 생성 등 모든 AI 호출이 이 모듈을 통합니다.

사용법:
    from app.services._resolve_ai import call_ai_json

    result = await call_ai_json(
        provider="gemini",           # "gemini" | "openai"
        model="gemini-2.5-flash",
        api_key="AIza...",
        system_prompt="You are...",
        user_content="Classify...",
    )
"""

import asyncio
import json
import logging
import random
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

MAX_RETRIES = 3
INITIAL_BACKOFF_S = 2.0
MAX_BACKOFF_S = 60.0
JITTER_FACTOR = 0.25

_RATE_LIMIT_INDICATORS = (
    "429",
    "resource_exhausted",
    "rate limit",
    "rate_limit",
    "quota",
    "too many requests",
)


def _is_rate_limit_error(exc: Exception) -> bool:
    msg = str(exc).lower()
    return any(indicator in msg for indicator in _RATE_LIMIT_INDICATORS)


async def call_ai_json(
    provider: str,
    model: str,
    api_key: str,
    system_prompt: str,
    user_content: str,
    temperature: float = 0.1,
    max_tokens: int = 1000,
) -> Optional[Dict[str, Any]]:
    """
    JSON 응답을 반환하는 AI 호출 (Gemini / OpenAI).
    Rate limit (429) 시 exponential backoff + jitter로 자동 재시도합니다.

    Returns:
        파싱된 JSON dict, 실패 시 None
    """
    provider = (provider or "gemini").lower().strip()

    if provider == "gemini":
        call_fn = _call_gemini
    elif provider == "openai":
        call_fn = _call_openai
    else:
        logger.error(f"Unknown AI provider: {provider}. Use 'gemini' or 'openai'.")
        return None

    last_exc: Optional[Exception] = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return await call_fn(model, api_key, system_prompt, user_content, temperature, max_tokens)
        except _RateLimitError as e:
            last_exc = e.original
            backoff = min(INITIAL_BACKOFF_S * (2 ** (attempt - 1)), MAX_BACKOFF_S)
            jitter = backoff * JITTER_FACTOR * random.random()
            wait = backoff + jitter
            logger.warning(
                f"Rate limited by {provider} (attempt {attempt}/{MAX_RETRIES}), "
                f"retrying in {wait:.1f}s: {e.original}"
            )
            await asyncio.sleep(wait)
        except _NonRetryableError as e:
            logger.error(f"{provider} call failed (model={model}): {e.original}")
            return None

    logger.error(
        f"{provider} call failed after {MAX_RETRIES} retries (model={model}): {last_exc}"
    )
    return None


class _RateLimitError(Exception):
    """Rate limit 에러를 래핑하여 재시도 로직에서 구분."""
    def __init__(self, original: Exception):
        self.original = original

class _NonRetryableError(Exception):
    """재시도 불필요한 에러 래핑."""
    def __init__(self, original: Exception):
        self.original = original


def _classify_and_raise(exc: Exception) -> None:
    if _is_rate_limit_error(exc):
        raise _RateLimitError(exc)
    raise _NonRetryableError(exc)


# ──────────────────────────────────────────────────────────────
# Gemini
# ──────────────────────────────────────────────────────────────

async def _call_gemini(
    model: str,
    api_key: str,
    system_prompt: str,
    user_content: str,
    temperature: float,
    max_tokens: int,
) -> Optional[Dict[str, Any]]:
    """
    google-genai SDK (v1+) 를 사용한 Gemini 호출.
    구 google-generativeai (deprecated 2025-11) 대신 google-genai 사용.
    """
    try:
        from google import genai
        from google.genai import types

        client = genai.Client(api_key=api_key)

        response = await client.aio.models.generate_content(
            model=model,
            contents=user_content,
            config=types.GenerateContentConfig(
                system_instruction=system_prompt,
                temperature=temperature,
                max_output_tokens=max_tokens,
                response_mime_type="application/json",
            ),
        )

        text = (response.text or "").strip()

        if text.startswith("```"):
            lines = text.splitlines()
            text = "\n".join(
                line for line in lines if not line.startswith("```")
            ).strip()

        return json.loads(text)

    except Exception as e:
        _classify_and_raise(e)


# ──────────────────────────────────────────────────────────────
# OpenAI
# ──────────────────────────────────────────────────────────────

async def _call_openai(
    model: str,
    api_key: str,
    system_prompt: str,
    user_content: str,
    temperature: float,
    max_tokens: int,
) -> Optional[Dict[str, Any]]:
    try:
        from openai import AsyncOpenAI

        client = AsyncOpenAI(api_key=api_key)
        response = await client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ],
            temperature=temperature,
            max_tokens=max_tokens,
            response_format={"type": "json_object"},
        )
        content = response.choices[0].message.content
        return json.loads(content)

    except Exception as e:
        _classify_and_raise(e)
