"""
Service Daemon - FastAPI Application

여러 유전체 분석 서비스를 통합 관리하는 범용 데몬입니다.
기존 nipt-daemon의 구조를 일반화하여 플러그인 기반으로 동작합니다.
"""

import os
import asyncio
import logging
from datetime import datetime
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Query
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse

from .config import settings
from .logging_config import setup_logging, setup_middleware
from .models import (
    OrderSubmitRequest, OrderSubmitResponse, OrderStatusResponse,
    OrderStatus, Job, QueueSummary
)
from .queue_manager import get_queue_manager
from .runner import get_runner
from .platform_client import get_platform_client
from .services import load_plugins, get_plugin, list_service_codes, get_all_plugins

logger = logging.getLogger(__name__)


# ─── Application Lifespan ──────────────────────────────────

@asynccontextmanager
async def lifespan(app: FastAPI):
    """애플리케이션 시작/종료 라이프사이클"""
    # Startup
    setup_logging()
    logger.info("=" * 60)
    logger.info(f"  Service Daemon Starting")
    logger.info(f"  Environment: {settings.app_env}")
    logger.info(f"  Max Concurrent Jobs: {settings.max_concurrent_jobs}")
    logger.info(f"  Enabled Services: {settings.enabled_service_list}")
    logger.info("=" * 60)

    # 서비스 플러그인 로딩
    load_plugins(settings.enabled_service_list)
    logger.info(f"Loaded services: {list_service_codes()}")

    # 워커 시작 (백그라운드)
    runner = get_runner()
    runner_task = asyncio.create_task(runner.start())

    yield

    # Shutdown
    logger.info("Shutting down Service Daemon...")
    runner._shutdown_event.set()
    try:
        await asyncio.wait_for(runner_task, timeout=30)
    except asyncio.TimeoutError:
        logger.warning("Runner shutdown timed out")
    logger.info("Service Daemon stopped")


# ─── FastAPI App ───────────────────────────────────────────

app = FastAPI(
    title="Service Daemon",
    description="Multi-service genomics pipeline daemon",
    version="1.0.0",
    lifespan=lifespan
)

# 미들웨어 설정
setup_middleware(app)

# 정적 파일
static_dir = os.path.join(os.path.dirname(__file__), "static")
if os.path.exists(static_dir):
    app.mount("/static", StaticFiles(directory=static_dir), name="static")


# ─── Health Check ──────────────────────────────────────────

@app.get("/health")
async def health():
    """헬스 체크"""
    return {
        "status": "healthy",
        "service": "service-daemon",
        "environment": settings.app_env,
        "timestamp": datetime.now().isoformat(),
        "registered_services": list_service_codes()
    }


# ─── Order Endpoints ──────────────────────────────────────

@app.post("/order/{service_code}/submit", response_model=OrderSubmitResponse)
async def submit_order(service_code: str, request: OrderSubmitRequest):
    """
    주문 접수 (서비스별).
    
    Platform에서 분석 주문을 접수합니다.
    service_code에 해당하는 플러그인이 등록되어 있어야 합니다.
    """
    # 서비스 플러그인 확인
    plugin = get_plugin(service_code)
    if not plugin:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown service_code: {service_code}. "
                   f"Available: {list_service_codes()}"
        )

    # 서비스 코드 일관성 확인
    if request.service_code != service_code:
        request.service_code = service_code

    # 파라미터 유효성 검사
    is_valid, error_msg = plugin.validate_params(request.params or {})
    if not is_valid:
        raise HTTPException(status_code=400, detail=f"Invalid params: {error_msg}")

    # 작업 디렉토리 결정
    work_dir = request.work_dir or datetime.now().strftime("%y%m%d")

    # Job 생성
    job = Job(
        order_id=request.order_id,
        service_code=service_code,
        sample_name=request.sample_name,
        work_dir=work_dir,
        fastq_r1_url=request.fastq_r1_url,
        fastq_r2_url=request.fastq_r2_url,
        fastq_r1_path=request.fastq_r1_path,
        fastq_r2_path=request.fastq_r2_path,
        params=request.params or {},
        priority=request.priority,
        callback_url=request.callback_url,
        # 경로 설정
        fastq_dir=os.path.join(settings.fastq_base_dir, work_dir, request.sample_name),
        analysis_dir=os.path.join(
            settings.analysis_base_dir, service_code, work_dir, request.sample_name
        ),
        output_dir=os.path.join(
            settings.output_base_dir, service_code, work_dir, request.sample_name
        ),
        log_dir=os.path.join(
            settings.log_base_dir, service_code, work_dir, request.sample_name
        )
    )

    # 큐에 추가
    queue_manager = get_queue_manager()
    queue_position = await queue_manager.enqueue(job)

    # 플랫폼에 접수 확인
    platform_client = get_platform_client()
    await platform_client.update_order_status(
        job.order_id, service_code,
        OrderStatus.QUEUED.value, 0,
        f"Order accepted, queue position: {queue_position}"
    )

    return OrderSubmitResponse(
        status="accepted",
        order_id=request.order_id,
        service_code=service_code,
        message=f"Order submitted to {plugin.display_name}",
        queue_position=queue_position
    )


@app.get("/order/{order_id}/status", response_model=OrderStatusResponse)
async def get_order_status(order_id: str):
    """주문 상태 조회"""
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    return OrderStatusResponse(
        order_id=job.order_id,
        service_code=job.service_code,
        status=job.status,
        progress=job.progress,
        message=job.message,
        created_at=job.created_at,
        updated_at=getattr(job, 'updated_at', None)
    )


@app.post("/order/{order_id}/cancel")
async def cancel_order(order_id: str):
    """주문 취소"""
    runner = get_runner()
    success = await runner.cancel_job(order_id)

    if success:
        return {"status": "cancelled", "order_id": order_id}
    else:
        raise HTTPException(
            status_code=404,
            detail=f"Cannot cancel order {order_id}: not found or not running"
        )


# ─── Queue Endpoints ──────────────────────────────────────

@app.get("/queue/summary", response_model=QueueSummary)
async def get_queue_summary():
    """큐 상태 요약"""
    queue_manager = get_queue_manager()
    return queue_manager.get_summary()


@app.get("/queue/status")
async def get_queue_status(
    service_code: str = Query(default=None, description="특정 서비스만 조회")
):
    """큐 상태 상세"""
    queue_manager = get_queue_manager()
    summary = queue_manager.get_summary()
    running_jobs = queue_manager.get_running_jobs(service_code)

    return {
        "summary": summary.model_dump(),
        "running_jobs": [
            {
                "order_id": j.order_id,
                "service_code": j.service_code,
                "sample_name": j.sample_name,
                "status": j.status.value,
                "progress": j.progress,
                "message": j.message,
                "started_at": j.started_at
            }
            for j in running_jobs
        ],
        "available_slots": queue_manager.available_slots,
        "max_concurrent": queue_manager.max_concurrent
    }


# ─── Service Info ──────────────────────────────────────────

@app.get("/services")
async def list_services():
    """등록된 서비스 목록"""
    plugins = get_all_plugins()
    return {
        "services": [
            {
                "service_code": code,
                "display_name": plugin.display_name,
                "progress_stages": plugin.get_progress_stages()
            }
            for code, plugin in plugins.items()
        ],
        "total": len(plugins)
    }


# ─── Dashboard ─────────────────────────────────────────────

@app.get("/")
async def dashboard():
    """간단한 대시보드 정보"""
    queue_manager = get_queue_manager()
    summary = queue_manager.get_summary()
    plugins = get_all_plugins()

    return {
        "service": "service-daemon",
        "version": "1.0.0",
        "environment": settings.app_env,
        "registered_services": [
            {"code": code, "name": plugin.display_name}
            for code, plugin in plugins.items()
        ],
        "queue": summary.model_dump(),
        "timestamp": datetime.now().isoformat()
    }
