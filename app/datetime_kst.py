"""KST (Asia/Seoul) timestamps for API responses, persisted jobs, and logs."""

from __future__ import annotations

from datetime import datetime
from zoneinfo import ZoneInfo

KST = ZoneInfo("Asia/Seoul")


def now_kst_iso() -> str:
    """ISO-8601 with +09:00 offset for JSON/job fields."""
    return datetime.now(KST).isoformat(timespec="seconds")


def now_kst_date_compact() -> str:
    """YYMMDD suffix for work_dir folders."""
    return datetime.now(KST).strftime("%y%m%d")


def now_kst_date_iso() -> str:
    """YYYY-MM-DD calendar date in Korea (e.g. report cover)."""
    return datetime.now(KST).strftime("%Y-%m-%d")
