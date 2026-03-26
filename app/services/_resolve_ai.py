"""
AI 제공자 추상화 레이어

Gemini / OpenAI 를 동일한 인터페이스로 호출합니다.
acmg.py, 미래의 report 생성 등 모든 AI 호출이 이 모듈을 통합니다.

사용법:
    from app.services._resolve_ai import call_ai_json

    result = await call_ai_json(
        provider="gemini",           # "gemini" | "openai"
        model="gemini-2.0-flash",
        api_key="AIza...",
        system_prompt="You are...",
        user_content="Classify...",
    )
"""

import json
import logging
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


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

    Returns:
        파싱된 JSON dict, 실패 시 None
    """
    provider = (provider or "gemini").lower().strip()

    if provider == "gemini":
        return await _call_gemini(model, api_key, system_prompt, user_content, temperature, max_tokens)
    elif provider == "openai":
        return await _call_openai(model, api_key, system_prompt, user_content, temperature, max_tokens)
    else:
        logger.error(f"Unknown AI provider: {provider}. Use 'gemini' or 'openai'.")
        return None


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

        # JSON 블록 추출 (마크다운 코드블록으로 감싸진 경우 대비)
        if text.startswith("```"):
            lines = text.splitlines()
            text = "\n".join(
                line for line in lines if not line.startswith("```")
            ).strip()

        return json.loads(text)

    except Exception as e:
        logger.error(f"Gemini call failed (model={model}): {e}")
        return None


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
        logger.error(f"OpenAI call failed (model={model}): {e}")
        return None
