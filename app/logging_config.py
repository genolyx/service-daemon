"""
Logging Configuration

기존 nipt-daemon의 logging_config.py를 일반화합니다.
"""

import logging
import time
from datetime import datetime
from colorlog import ColoredFormatter
from fastapi import Request

from .config import settings
from .datetime_kst import KST

log = logging.getLogger(__name__)


def _kst_log_converter(seconds: float):
    return datetime.fromtimestamp(seconds, tz=KST).timetuple()


def setup_logging():
    """로깅 설정 초기화"""
    logging_level = getattr(logging, settings.log_level.upper(), logging.INFO)

    formatter = ColoredFormatter(
        "%(log_color)s%(asctime)s %(levelname)s %(filename)s:%(lineno)d ─ %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        reset=True,
        log_colors={
            "DEBUG": "cyan",
            "INFO": "white",
            "WARNING": "yellow",
            "ERROR": "red",
            "CRITICAL": "red",
        },
        secondary_log_colors={},
        style="%",
    )
    formatter.converter = _kst_log_converter

    logging.basicConfig(
        level=logging_level,
        format="%(asctime)s %(levelname)s %(name)s [%(filename)s:%(lineno)d] %(message)s"
    )

    handler = logging.StreamHandler()
    handler.setLevel(logging_level)
    handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.setLevel(logging_level)
    root_logger.handlers = [handler]

    # httpx 로그 레벨 설정
    logging.getLogger("httpx").setLevel(logging.WARNING)


def setup_middleware(app):
    """HTTP 요청 로깅 미들웨어 (한 줄 요약)"""

    @app.middleware("http")
    async def log_requests(request: Request, call_next):
        skip_logging_paths = [
            "/health",
            "/queue/summary",
            "/queue/status",
            "/favicon.ico"
        ]

        skip_user_agents = [
            "Hello World", "curl", "wget", "bot", "crawler", "spider", "scanner"
        ]

        should_skip_path = any(request.url.path == path for path in skip_logging_paths)
        user_agent = request.headers.get("user-agent", "").lower()
        should_skip_ua = any(pattern.lower() in user_agent for pattern in skip_user_agents)
        is_external_root = (
            request.url.path == "/" and
            request.client.host not in ["127.0.0.1", "localhost", "::1"]
        )

        should_skip = should_skip_path or should_skip_ua or is_external_root

        start_time = time.time()
        response = await call_next(request)
        process_time = time.time() - start_time

        if not should_skip:
            path_qs = request.url.path
            if request.url.query:
                path_qs = f"{path_qs}?{request.url.query}"
            peer = (
                f"{request.client.host}:{request.client.port}"
                if request.client else "-"
            )
            log.info(
                "%s %s -> %s %.3fs %s",
                request.method,
                path_qs,
                response.status_code,
                process_time,
                peer,
            )
        else:
            log.debug(
                "skip access log: %s %s from %s",
                request.method,
                request.url.path,
                request.client.host if request.client else "-",
            )

        return response
