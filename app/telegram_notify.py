"""
Optional Telegram notifications for order lifecycle (runner hooks).

Uses Bot API sendMessage (outbound HTTPS). No inbound webhook required.
"""

from __future__ import annotations

import asyncio
import logging
from typing import List

import httpx

from .config import settings
from .models import Job

logger = logging.getLogger(__name__)


def _parse_chat_ids(raw: str) -> List[str]:
    return [p.strip() for p in raw.replace(";", ",").split(",") if p.strip()]


async def send_order_telegram(event: str, job: Job) -> None:
    """
    event: started | completed | failed | cancelled
    """
    if not settings.telegram_notify_enabled:
        return
    token = (settings.telegram_bot_token or "").strip()
    raw_ids = (settings.telegram_chat_ids or "").strip()
    if not token or not raw_ids:
        return

    label_map = {
        "started": "Started",
        "completed": "Completed",
        "failed": "Fail",
        "cancelled": "Cancelled",
    }
    label = label_map.get(event)
    if not label:
        logger.warning("Unknown telegram event: %s", event)
        return

    text = f"[{job.order_id}] {label}"
    url = f"https://api.telegram.org/bot{token}/sendMessage"
    chat_ids = _parse_chat_ids(raw_ids)

    try:
        async with httpx.AsyncClient(timeout=15.0) as client:
            for cid in chat_ids:
                r = await client.post(
                    url,
                    json={
                        "chat_id": cid,
                        "text": text,
                        "disable_web_page_preview": True,
                    },
                )
                if r.status_code != 200:
                    logger.warning(
                        "Telegram sendMessage non-200 for chat_id=%s: %s %s",
                        cid,
                        r.status_code,
                        (r.text or "")[:500],
                    )
                else:
                    body = r.json()
                    if not body.get("ok"):
                        logger.warning(
                            "Telegram API ok=false for chat_id=%s: %s",
                            cid,
                            body,
                        )
    except Exception as e:
        logger.warning("Telegram notification failed: %s", e)


def schedule_order_telegram(event: str, job: Job) -> None:
    """Fire-and-forget so Telegram latency/outages never block the runner."""

    async def _run() -> None:
        try:
            await send_order_telegram(event, job)
        except Exception as e:
            logger.warning("Telegram task error: %s", e)

    try:
        asyncio.get_running_loop().create_task(_run())
    except RuntimeError:
        logger.debug("No running event loop; skipping Telegram schedule")
