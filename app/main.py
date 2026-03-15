"""
Service Daemon - FastAPI Application

여러 유전체 분석 서비스를 통합 관리하는 범용 데몬입니다.
기존 nipt-daemon의 구조를 일반화하여 플러그인 기반으로 동작합니다.
"""

import os
import json
import asyncio
import glob
import logging
from datetime import datetime
from contextlib import asynccontextmanager
from typing import Dict, Any, List, Optional

from fastapi import FastAPI, HTTPException, Query, Body
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware

from .config import settings
from .logging_config import setup_logging, setup_middleware
from .models import (
    OrderSubmitRequest, OrderSubmitResponse, OrderStatusResponse,
    OrderStatus, Job, QueueSummary, OutputFile,
    ReportGenerateRequest, ReportGenerateResponse,
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
    version="2.0.0",
    lifespan=lifespan
)

# CORS 미들웨어 (테스트 Portal에서 로컬 API 호출 허용)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
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

    # 디렉토리 경로 설정 (서비스별 구조)
    base = os.path.join(settings.base_dir, "carrier-screening") \
        if service_code == "carrier_screening" \
        else settings.base_dir

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
        fastq_dir=os.path.join(base, "fastq", work_dir, request.sample_name),
        analysis_dir=os.path.join(base, "analysis", work_dir, request.sample_name),
        output_dir=os.path.join(base, "output", work_dir, request.sample_name),
        log_dir=os.path.join(base, "log", work_dir, request.sample_name),
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


# ─── Report Generation Endpoint ───────────────────────────

@app.post("/order/{order_id}/report", response_model=ReportGenerateResponse)
async def generate_report(order_id: str, request: ReportGenerateRequest):
    """
    리뷰어 확정 후 최종 리포트를 생성합니다.

    Portal에서 리뷰어가 변이를 확정하고 코멘트를 작성한 후,
    이 엔드포인트를 호출하여 report.json + 다국어 PDF를 생성합니다.
    생성된 리포트는 자동으로 Platform에 업로드됩니다.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    plugin = get_plugin(job.service_code)
    if not plugin:
        raise HTTPException(
            status_code=400,
            detail=f"No plugin for service: {job.service_code}"
        )

    # generate_report 메서드가 있는지 확인
    if not hasattr(plugin, 'generate_report'):
        raise HTTPException(
            status_code=400,
            detail=f"Service {job.service_code} does not support report generation"
        )

    # 리포트 생성 (patient_info, partner_info, languages 전달)
    success = await plugin.generate_report(
        job=job,
        confirmed_variants=request.confirmed_variants,
        reviewer_info=request.reviewer_info,
        patient_info=request.patient_info,
        partner_info=request.partner_info,
        languages=request.languages,
    )

    if not success:
        raise HTTPException(
            status_code=500,
            detail="Report generation failed"
        )

    # 생성된 리포트 파일을 Platform에 업로드
    output_dir = job.output_dir
    report_files = []

    # report.json
    report_json = os.path.join(output_dir, "report.json")
    if os.path.exists(report_json):
        report_files.append(OutputFile(
            file_path=report_json,
            file_type="report_json",
            file_name="report.json",
            content_type="application/json",
        ))

    # 다국어 PDF 파일 (Report_*.pdf 패턴)
    pdf_files = glob.glob(os.path.join(output_dir, "Report_*.pdf"))
    for pdf_path in pdf_files:
        report_files.append(OutputFile(
            file_path=pdf_path,
            file_type="report_pdf",
            file_name=os.path.basename(pdf_path),
            content_type="application/pdf",
        ))

    # 단일 report.pdf (내장 템플릿 사용 시)
    single_pdf = os.path.join(output_dir, "report.pdf")
    if os.path.exists(single_pdf) and single_pdf not in [f.file_path for f in report_files]:
        report_files.append(OutputFile(
            file_path=single_pdf,
            file_type="report_pdf",
            file_name="report.pdf",
            content_type="application/pdf",
        ))

    # HTML 리포트 파일
    html_files = glob.glob(os.path.join(output_dir, "Report_*.html"))
    for html_path in html_files:
        report_files.append(OutputFile(
            file_path=html_path,
            file_type="report_html",
            file_name=os.path.basename(html_path),
            content_type="text/html",
        ))

    # Platform에 업로드
    platform_client = get_platform_client()
    uploaded_count = 0
    if report_files:
        upload_results = await platform_client.upload_all_outputs(
            order_id, job.service_code, report_files
        )
        uploaded_count = sum(
            1 for r in upload_results.values()
            if r.status.value == "SUCCESS"
        )

    logger.info(
        f"Report generation complete for {order_id}: "
        f"{len(report_files)} files generated, {uploaded_count} uploaded"
    )

    return ReportGenerateResponse(
        status="success",
        order_id=order_id,
        service_code=job.service_code,
        report_files=[os.path.basename(f.file_path) for f in report_files],
        uploaded_count=uploaded_count,
        message=f"Report generated: {len(report_files)} file(s), {uploaded_count} uploaded",
    )


# ─── Order List & Result Endpoints ────────────────────────

@app.get("/orders")
async def list_orders(
    service_code: str = Query(default=None, description="특정 서비스만 조회"),
    status: str = Query(default=None, description="특정 상태만 조회"),
):
    """
    전체 주문 목록을 조회합니다.

    Portal의 Track Orders 뷰에서 사용합니다.
    """
    queue_manager = get_queue_manager()
    all_jobs = []

    # 모든 Job 소스에서 수집
    for order_id, job in queue_manager._jobs.items():
        all_jobs.append(job)
    for order_id, job in queue_manager._running_jobs.items():
        if order_id not in {j.order_id for j in all_jobs}:
            all_jobs.append(job)
    for order_id, job in queue_manager._completed_jobs.items():
        if order_id not in {j.order_id for j in all_jobs}:
            all_jobs.append(job)

    # 필터링
    if service_code:
        all_jobs = [j for j in all_jobs if j.service_code == service_code]
    if status:
        all_jobs = [j for j in all_jobs if j.status.value == status]

    # 최신순 정렬
    all_jobs.sort(key=lambda j: j.created_at or "", reverse=True)

    return {
        "orders": [
            {
                "order_id": j.order_id,
                "service_code": j.service_code,
                "sample_name": j.sample_name,
                "status": j.status.value,
                "progress": j.progress,
                "message": j.message,
                "created_at": j.created_at,
                "started_at": j.started_at,
                "completed_at": j.completed_at,
                "priority": j.priority,
            }
            for j in all_jobs
        ],
        "total": len(all_jobs),
    }


@app.get("/order/{order_id}/result")
async def get_order_result(order_id: str):
    """
    주문의 result.json을 조회합니다.

    파이프라인 완료 후 process_results에서 생성된 result.json을 반환합니다.
    Portal의 Variant Review 페이지에서 사용합니다.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    output_dir = job.output_dir
    if not output_dir:
        raise HTTPException(
            status_code=404,
            detail=f"Output directory not set for order: {order_id}"
        )

    result_json_path = os.path.join(output_dir, "result.json")
    if not os.path.exists(result_json_path):
        raise HTTPException(
            status_code=404,
            detail=f"result.json not found for order: {order_id}. "
                   f"Analysis may not be complete yet."
        )

    with open(result_json_path, "r", encoding="utf-8") as f:
        result_data = json.load(f)

    return result_data


@app.get("/order/{order_id}/files")
async def list_order_files(order_id: str):
    """
    주문의 출력 파일 목록을 조회합니다.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    output_dir = job.output_dir
    if not output_dir or not os.path.isdir(output_dir):
        return {"files": [], "total": 0}

    files = []
    for fname in os.listdir(output_dir):
        fpath = os.path.join(output_dir, fname)
        if os.path.isfile(fpath):
            files.append({
                "name": fname,
                "size": os.path.getsize(fpath),
                "type": _guess_file_type(fname),
            })

    return {"files": files, "total": len(files)}


@app.get("/order/{order_id}/file/{filename}")
async def download_order_file(order_id: str, filename: str):
    """
    주문의 특정 출력 파일을 다운로드합니다.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    output_dir = job.output_dir
    if not output_dir:
        raise HTTPException(status_code=404, detail="Output directory not set")

    file_path = os.path.join(output_dir, filename)
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail=f"File not found: {filename}")

    return FileResponse(
        path=file_path,
        filename=filename,
        media_type=_guess_content_type(filename),
    )


def _guess_file_type(filename: str) -> str:
    """파일명으로 파일 유형을 추정합니다."""
    ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""
    type_map = {
        "json": "json", "tsv": "tsv", "csv": "csv",
        "pdf": "pdf", "html": "html",
        "png": "image", "jpg": "image", "jpeg": "image",
        "vcf": "vcf", "bed": "bed",
    }
    return type_map.get(ext, "other")


def _guess_content_type(filename: str) -> str:
    """파일명으로 Content-Type을 추정합니다."""
    ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""
    ct_map = {
        "json": "application/json",
        "tsv": "text/tab-separated-values",
        "csv": "text/csv",
        "pdf": "application/pdf",
        "html": "text/html",
        "png": "image/png",
        "jpg": "image/jpeg",
        "jpeg": "image/jpeg",
    }
    return ct_map.get(ext, "application/octet-stream")


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
        "version": "2.0.0",
        "environment": settings.app_env,
        "registered_services": [
            {"code": code, "name": plugin.display_name}
            for code, plugin in plugins.items()
        ],
        "queue": summary.model_dump(),
        "timestamp": datetime.now().isoformat()
    }


# ─── Test/Mock Endpoints (Development Only) ──────────────

@app.post("/test/inject-mock-job")
async def inject_mock_job(
    order_id: str = Body(default="CS-TEST-001"),
    sample_name: str = Body(default="SAMPLE_001"),
    output_dir: str = Body(default="/tmp/carrier-screening/output/test/SAMPLE_001"),
    status: str = Body(default="COMPLETED"),
):
    """
    [DEV ONLY] 테스트용 mock Job을 QueueManager에 주입합니다.
    실제 파이프라인 없이 Portal의 Review/Report 기능을 테스트할 수 있습니다.
    """
    if settings.app_env not in ("DEV", "TEST", "PROD"):
        raise HTTPException(status_code=403, detail="Mock injection only available in DEV/TEST")

    queue_manager = get_queue_manager()

    # 이미 존재하면 삭제
    existing = queue_manager.get_job(order_id)
    if existing:
        return {"status": "exists", "order_id": order_id, "message": "Job already exists"}

    job = Job(
        order_id=order_id,
        service_code="carrier_screening",
        sample_name=sample_name,
        work_dir="test",
        fastq_r1_path="/tmp/test_R1.fastq.gz",
        fastq_r2_path="/tmp/test_R2.fastq.gz",
        params={},
        priority="normal",
        fastq_dir="/tmp/carrier-screening/fastq/test/" + sample_name,
        analysis_dir="/tmp/carrier-screening/analysis/test/" + sample_name,
        output_dir=output_dir,
        log_dir="/tmp/carrier-screening/log/test/" + sample_name,
    )

    # 상태 설정
    status_enum = OrderStatus(status) if status in [s.value for s in OrderStatus] else OrderStatus.COMPLETED
    job.status = status_enum
    job.created_at = datetime.now().isoformat()
    job.progress = 100 if status_enum == OrderStatus.COMPLETED else 0

    if status_enum == OrderStatus.COMPLETED:
        job.completed_at = datetime.now().isoformat()
        job.started_at = datetime.now().isoformat()
        job.message = "Mock job - analysis complete"
        queue_manager._completed_jobs[order_id] = job
    elif status_enum == OrderStatus.RUNNING:
        job.started_at = datetime.now().isoformat()
        job.message = "Mock job - running"
        queue_manager._running_jobs[order_id] = job
    else:
        job.message = "Mock job - queued"
        queue_manager._jobs[order_id] = job

    logger.info(f"[TEST] Injected mock job: {order_id} ({status})")

    return {
        "status": "injected",
        "order_id": order_id,
        "sample_name": sample_name,
        "job_status": job.status.value,
        "output_dir": output_dir,
    }
