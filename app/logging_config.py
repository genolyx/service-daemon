"""
Logging Configuration

기존 nipt-daemon의 logging_config.py를 일반화합니다.
"""

import logging
import time
import json
from datetime import datetime
from colorlog import ColoredFormatter
from fastapi import Request

from .config import settings


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
    """HTTP 요청 로깅 미들웨어"""

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

        if not should_skip:
            await pretty_print_request(request)

        response = await call_next(request)
        process_time = time.time() - start_time

        if not should_skip:
            print(f"  Response Time: {process_time:.3f}s | Status: {response.status_code}")
            print("=" * 80 + "\n")
        else:
            logger = logging.getLogger(__name__)
            logger.debug(
                f"Skipped logging: {request.method} {request.url.path} "
                f"from {request.client.host}"
            )

        return response


async def pretty_print_request(request: Request):
    """Request 정보 출력"""
    print("\n" + "=" * 80)
    print(f" REQUEST LOG - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f" Client: {request.client.host}:{request.client.port}")
    print(f" {request.method} {request.url}")
    print("-" * 80)

    important_headers = ['authorization', 'content-type', 'user-agent', 'accept']
    print(" HEADERS:")
    for name, value in request.headers.items():
        if name.lower() in important_headers:
            if name.lower() == 'authorization':
                print(f"   {name}: Bearer ***masked***")
            else:
                print(f"   {name}: {value}")

    try:
        body = await request.body()
        if body:
            print("\n BODY:")
            try:
                json_body = json.loads(body.decode('utf-8'))
                print(json.dumps(json_body, indent=3, ensure_ascii=False))
            except json.JSONDecodeError:
                body_str = body.decode('utf-8')
                if len(body_str) > 200:
                    print(f"   {body_str[:200]}... (truncated)")
                else:
                    print(f"   {body_str}")
        else:
            print("\n BODY: (empty)")
    except Exception as e:
        print(f"\n BODY: Error reading - {e}")
