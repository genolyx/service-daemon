"""
Service Daemon - FastAPI Application

여러 유전체 분석 서비스를 통합 관리하는 범용 데몬입니다.
기존 nipt-daemon의 구조를 일반화하여 플러그인 기반으로 동작합니다.
"""

import os
import json
import asyncio
import glob
import hmac
import logging
from contextlib import asynccontextmanager
from typing import Dict, Any, List, Optional

from fastapi import FastAPI, HTTPException, Query, Body, Request
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse, FileResponse, PlainTextResponse, RedirectResponse
from fastapi.middleware.cors import CORSMiddleware

from .config import settings
from .datetime_kst import now_kst_iso, now_kst_date_compact
from .logging_config import setup_logging, setup_middleware
from .models import (
    OrderSubmitRequest, OrderSubmitResponse, OrderSaveResponse, OrderStatusResponse,
    OrderUpdateRequest, OrderUpdateResponse,
    OrderStatus, Job, QueueSummary, OutputFile,
    ReportGenerateRequest, ReportGenerateResponse,
    GeneKnowledgeSaveRequest,
    VariantKnowledgeSaveRequest,
    UpdateFastqPathsRequest,
)
from .queue_manager import get_queue_manager
from .order_store import ingest_report_json_from_disk
from .annotation_resources import annotation_resource_report
from .runner import get_runner
from .platform_client import get_platform_client
from .services import load_plugins, get_plugin, list_service_codes, get_all_plugins

logger = logging.getLogger(__name__)

_FASTQ_NAME_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _fastq_root_real() -> str:
    return os.path.realpath(settings.fastq_base_dir)


def _fastq_root_for_service(service_code: Optional[str]) -> str:
    """포털 browse / FASTQ 경로 검증용 루트 (서비스별)."""
    sc = (service_code or "").strip().lower().replace("-", "_")
    if sc == "sgnipt":
        return os.path.realpath(settings.sgnipt_fastq_root)
    if sc == "carrier_screening":
        return os.path.realpath(settings.carrier_screening_fastq_dir)
    return _fastq_root_real()


def _normalize_rel_path(rel: str) -> str:
    rel = (rel or "").replace("\\", "/").strip().lstrip("/")
    parts: List[str] = []
    for p in rel.split("/"):
        if not p or p == ".":
            continue
        if p == "..":
            raise HTTPException(status_code=400, detail="Invalid path segment")
        parts.append(p)
    return "/".join(parts)


def _safe_join_under_fastq(rel: str, service_code: Optional[str] = None) -> str:
    """서비스별 FASTQ 베이스 아래 절대 경로. rel 빈 문자열이면 루트."""
    root = _fastq_root_for_service(service_code)
    norm = _normalize_rel_path(rel)
    target = os.path.realpath(os.path.join(root, norm) if norm else root)
    root_sep = root if root.endswith(os.sep) else root + os.sep
    if target != root and not target.startswith(root_sep):
        raise HTTPException(status_code=400, detail="Path escapes FASTQ root")
    return target


def _is_fastq_filename(name: str) -> bool:
    lower = name.lower()
    return any(lower.endswith(s) for s in _FASTQ_NAME_SUFFIXES)


def _validate_optional_fastq_file(abs_path: str, service_code: str) -> str:
    full = os.path.realpath(abs_path.strip())
    if not os.path.isfile(full):
        raise HTTPException(status_code=400, detail=f"Not a file: {abs_path}")
    root = _fastq_root_for_service(service_code)
    root_sep = root if root.endswith(os.sep) else root + os.sep
    if full != root and not full.startswith(root_sep):
        raise HTTPException(
            status_code=400,
            detail=f"File must be under this service FASTQ directory: {root}",
        )
    return full


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
    logger.info(f"  Orders DB (SQLite): {settings.resolved_orders_db_path}")
    logger.info("=" * 60)

    # 서비스 플러그인 로딩
    load_plugins(settings.enabled_service_list)
    logger.info(f"Loaded services: {list_service_codes()}")

    # Literature 캐시 DB 초기화 (PubMed 검색 결과 영구 캐시)
    if settings.literature_enabled:
        try:
            from .services.carrier_screening.literature import init_literature_db
            lit_db = settings.resolved_literature_db_path
            os.makedirs(os.path.dirname(lit_db), exist_ok=True)
            init_literature_db(lit_db)
            logger.info(f"Literature DB (SQLite): {lit_db}")
        except Exception as e:
            logger.warning(f"Literature DB init failed (non-fatal): {e}")

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


def _api_key_configured() -> bool:
    k = settings.api_key
    return bool(k and str(k).strip())


def _header_matches_api_key(header_value: str, key: str) -> bool:
    if not header_value or not key:
        return False
    if len(header_value) != len(key):
        return False
    return hmac.compare_digest(header_value, key)


@app.middleware("http")
async def api_access_key_guard(request: Request, call_next):
    """
    선택적 인바운드 API 키 (추후 외부 Portal / 자동화 클라이언트 대비).
    API_KEY 미설정 시 통과. 설정 시 portal·헬스·정적·OPTIONS 는 예외.
    허용: Authorization: Bearer <api_key> 또는 X-API-Key: <api_key>
    """
    if not _api_key_configured():
        return await call_next(request)

    if request.method == "OPTIONS":
        return await call_next(request)

    path = request.url.path
    if path == "/health" or path.startswith("/portal") or path.startswith("/static"):
        return await call_next(request)

    key = str(settings.api_key).strip()
    auth = request.headers.get("Authorization") or ""
    expected_bearer = f"Bearer {key}"
    ok = False
    if len(auth) == len(expected_bearer):
        ok = hmac.compare_digest(auth, expected_bearer)
    if not ok:
        xk = request.headers.get("X-API-Key") or ""
        ok = _header_matches_api_key(xk, key)

    if ok:
        return await call_next(request)

    logger.warning("Forbidden request: invalid or missing API key (Authorization / X-API-Key)")
    return JSONResponse(status_code=403, content={"detail": "Forbidden"})


# 정적 파일
static_dir = os.path.join(os.path.dirname(__file__), "static")
if os.path.exists(static_dir):
    app.mount("/static", StaticFiles(directory=static_dir), name="static")

# Test Portal (carrier-screening, sgNIPT 연동 테스트용)
portal_dir = os.path.join(os.path.dirname(__file__), "..", "portal")
if os.path.exists(portal_dir):
    @app.get("/portal", include_in_schema=False)
    def portal_redirect():
        return RedirectResponse("/portal/")
    app.mount("/portal", StaticFiles(directory=portal_dir, html=True), name="portal")


# ─── Health Check ──────────────────────────────────────────

@app.get("/health")
async def health(
    annotation: bool = Query(
        default=False,
        description="true면 ClinVar/gnomAD/snpEff 등 주석용 리소스 존재 여부 요약 포함",
    ),
):
    """헬스 체크"""
    fq_default = _fastq_root_real()
    body: Dict[str, Any] = {
        "status": "healthy",
        "service": "service-daemon",
        "environment": settings.app_env,
        "timestamp": now_kst_iso(),
        "registered_services": list_service_codes(),
        "fastq_base_dir": fq_default,
        "fastq_root_exists": os.path.isdir(fq_default),
        "fastq_browse_roots": {
            "sgnipt": _fastq_root_for_service("sgnipt"),
            "carrier_screening": _fastq_root_for_service("carrier_screening"),
            "default": fq_default,
        },
    }
    if annotation:
        body["annotation_resources"] = annotation_resource_report()
    return body


# ─── Order Endpoints ──────────────────────────────────────

def _build_job_from_submit_request(service_code: str, request: OrderSubmitRequest) -> Job:
    """Submit/Save 공통 Job 생성."""
    work_dir = request.work_dir or now_kst_date_compact()
    if service_code == "sgnipt":
        root = settings.sgnipt_job_root
        oid = (request.order_id or "").strip()
        return Job(
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
            fastq_dir=os.path.join(root, "fastq", work_dir, oid),
            analysis_dir=os.path.join(root, "analysis", work_dir, oid),
            output_dir=os.path.join(root, "output", work_dir, oid),
            log_dir=os.path.join(root, "log", work_dir, oid),
        )
    if service_code == "carrier_screening":
        work_root = settings.carrier_screening_work_root
        return Job(
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
            fastq_dir=os.path.join(work_root, "fastq", work_dir, request.sample_name),
            analysis_dir=os.path.join(work_root, "analysis", work_dir, request.sample_name),
            output_dir=os.path.join(work_root, "output", work_dir, request.sample_name),
            log_dir=os.path.join(work_root, "log", work_dir, request.sample_name),
        )
    base = settings.base_dir
    return Job(
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
        fastq_dir=os.path.join(base, "fastq", work_dir, request.sample_name),
        analysis_dir=os.path.join(base, "analysis", work_dir, request.sample_name),
        output_dir=os.path.join(base, "output", work_dir, request.sample_name),
        log_dir=os.path.join(base, "log", work_dir, request.sample_name),
    )


def _merge_order_params(old: Dict[str, Any], patch: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(old or {})
    for k, v in patch.items():
        if k in ("carrier", "nipt") and isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = {**out[k], **v}
        else:
            out[k] = v
    return out


def _order_merge_submit(job: Job, patch: OrderUpdateRequest) -> OrderSubmitRequest:
    data = patch.model_dump(exclude_unset=True)
    new_oid = data.pop("order_id", None)
    params_patch = data.pop("params", None)
    merged_params = dict(job.params or {})
    if params_patch is not None:
        merged_params = _merge_order_params(merged_params, params_patch)
    base: Dict[str, Any] = {
        "order_id": job.order_id,
        "service_code": job.service_code,
        "sample_name": job.sample_name,
        "work_dir": job.work_dir,
        "fastq_r1_url": job.fastq_r1_url,
        "fastq_r2_url": job.fastq_r2_url,
        "fastq_r1_path": job.fastq_r1_path,
        "fastq_r2_path": job.fastq_r2_path,
        "params": merged_params,
        "priority": job.priority or "normal",
        "callback_url": job.callback_url,
    }
    for k, v in data.items():
        if v is not None:
            base[k] = v
    if new_oid is not None and str(new_oid).strip():
        base["order_id"] = str(new_oid).strip()
    return OrderSubmitRequest(**base)


@app.patch("/order/{order_id}", response_model=OrderUpdateResponse)
async def update_order(order_id: str, request: OrderUpdateRequest):
    """
    SAVED / FAILED / CANCELLED / COMPLETED / REPORT_READY 주문 필드 수정.
    (QUEUED·실행 중 등은 불가 — 완료 후 메타 수정·재실행 전에 변경 가능)
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    if job.status not in (
        OrderStatus.SAVED,
        OrderStatus.FAILED,
        OrderStatus.CANCELLED,
        OrderStatus.COMPLETED,
        OrderStatus.REPORT_READY,
    ):
        raise HTTPException(
            status_code=400,
            detail=(
                f"Order {order_id} is {job.status.value}; "
                "only SAVED, FAILED, CANCELLED, COMPLETED, or REPORT_READY can be edited"
            ),
        )
    submit = _order_merge_submit(job, request)
    plugin = get_plugin(submit.service_code)
    if not plugin:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown service_code: {submit.service_code}",
        )
    is_valid, error_msg = plugin.validate_params(submit.params or {})
    if not is_valid:
        raise HTTPException(status_code=400, detail=f"Invalid params: {error_msg}")
    new_job = _build_job_from_submit_request(submit.service_code, submit)
    new_job.status = job.status
    new_job.created_at = job.created_at
    new_job.updated_at = now_kst_iso()
    if job.status in (OrderStatus.FAILED, OrderStatus.CANCELLED):
        new_job.progress = 0
        new_job.message = ""
        new_job.error_log = None
        new_job.completed_at = None
        new_job.started_at = None
    elif job.status == OrderStatus.SAVED:
        new_job.progress = 0
        new_job.message = ""
    elif job.status in (OrderStatus.COMPLETED, OrderStatus.REPORT_READY):
        new_job.progress = job.progress
        new_job.message = job.message
        new_job.completed_at = job.completed_at
        new_job.started_at = job.started_at
        new_job.error_log = job.error_log
        new_job.pid = job.pid
        new_job.exit_code = job.exit_code
        new_job.duration = job.duration
    try:
        await queue_manager.replace_edited_job(new_job, previous_order_id=order_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return OrderUpdateResponse(
        status="updated",
        order_id=new_job.order_id,
        message="Order updated",
    )


@app.post("/order/{service_code}/save", response_model=OrderSaveResponse)
async def save_order(service_code: str, request: OrderSubmitRequest):
    """
    주문만 저장 (SAVED). 큐에 넣지 않음 — Portal에서 폼 유실 방지 후 More → Submit.
    """
    plugin = get_plugin(service_code)
    if not plugin:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown service_code: {service_code}. "
                   f"Available: {list_service_codes()}"
        )
    if request.service_code != service_code:
        request.service_code = service_code
    is_valid, error_msg = plugin.validate_params(request.params or {})
    if not is_valid:
        raise HTTPException(status_code=400, detail=f"Invalid params: {error_msg}")

    queue_manager = get_queue_manager()
    existing = queue_manager.get_job(request.order_id)
    if existing and existing.status != OrderStatus.SAVED:
        raise HTTPException(
            status_code=409,
            detail=f"Order {request.order_id} already exists with status {existing.status.value}",
        )

    job = _build_job_from_submit_request(service_code, request)
    await queue_manager.save_job(job)
    return OrderSaveResponse(
        status="saved",
        order_id=request.order_id,
        service_code=service_code,
        message=f"Order saved (not queued). Use POST /order/{{order_id}}/start to run.",
    )


@app.post("/order/{order_id}/start", response_model=OrderSubmitResponse)
async def start_saved_order(order_id: str):
    """SAVED 주문을 큐에 넣거나, FAILED/CANCELLED 주문을 재시도로 다시 큐에 넣습니다."""
    queue_manager = get_queue_manager()
    try:
        job, queue_position = await queue_manager.start_saved_job(order_id)
    except KeyError:
        raise HTTPException(
            status_code=404,
            detail=f"No startable order: {order_id} (need SAVED, FAILED, CANCELLED, COMPLETED, or REPORT_READY)",
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    platform_client = get_platform_client()
    try:
        await platform_client.update_order_status(
            job.order_id, job.service_code,
            OrderStatus.QUEUED.value, 0,
            f"Order accepted, queue position: {queue_position}"
        )
    except Exception as e:
        logger.warning(f"Platform status update after start failed: {e}")

    plugin = get_plugin(job.service_code)
    return OrderSubmitResponse(
        status="accepted",
        order_id=order_id,
        service_code=job.service_code,
        message=f"Order started ({plugin.display_name if plugin else job.service_code})",
        queue_position=queue_position,
    )


@app.post("/order/{order_id}/reprocess-results")
async def reprocess_order_results(order_id: str):
    """
    Carrier screening: Nextflow 없이 process_results 만 다시 실행합니다 (Portal «Reprocess only»).
    전체를 FASTQ부터 돌리려면 POST /order/{id}/start («Force Run»)를 사용하세요.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    if job.service_code != "carrier_screening":
        raise HTTPException(
            status_code=400,
            detail="reprocess-results is only supported for carrier_screening",
        )

    if job.status not in (
        OrderStatus.COMPLETED,
        OrderStatus.REPORT_READY,
        OrderStatus.FAILED,
    ):
        raise HTTPException(
            status_code=400,
            detail=(
                f"Order must be COMPLETED, REPORT_READY, or FAILED (current: {job.status.value}). "
                "Stop the run first if it is active; use POST /order/{id}/start to queue a full run."
            ),
        )

    plugin = get_plugin(job.service_code)
    if not plugin:
        raise HTTPException(status_code=400, detail="No plugin for carrier_screening")

    from .services.carrier_screening.layout_norm import apply_carrier_layout_directories

    if apply_carrier_layout_directories(job):
        await queue_manager.persist_job(job)

    completion_ok = await plugin.check_completion(job)
    if not completion_ok:
        raise HTTPException(
            status_code=400,
            detail="No VCF found for this order — cannot reprocess (check analysis/output paths).",
        )

    logger.info("[reprocess-results] Starting process_results for %s", order_id)
    process_ok = await plugin.process_results(job)
    if not process_ok:
        raise HTTPException(
            status_code=500,
            detail="process_results failed — see daemon logs",
        )

    await queue_manager.finalize_reprocess_results(job)
    logger.info("[reprocess-results] Done for %s", order_id)

    return {
        "status": "ok",
        "order_id": order_id,
        "message": "Reprocessed (annotation/QC/result.json updated)",
    }


@app.post("/order/{order_id}/stop")
async def stop_order(order_id: str):
    """
    SAVED 삭제 / QUEUED 취소 요청 / RUNNING 파이프라인 중단.
    """
    queue_manager = get_queue_manager()
    runner = get_runner()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    if job.status == OrderStatus.SAVED:
        ok = await queue_manager.delete_saved(order_id)
        if not ok:
            raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
        return {"status": "deleted", "order_id": order_id, "message": "Saved draft removed"}

    if job.status == OrderStatus.QUEUED:
        ok = await queue_manager.request_cancel_queued(order_id)
        if not ok:
            raise HTTPException(
                status_code=400,
                detail=f"Cannot cancel queued order {order_id}",
            )
        return {"status": "cancelling", "order_id": order_id, "message": "Will cancel when dequeued"}

    if job.status in (
        OrderStatus.RUNNING,
        OrderStatus.DOWNLOADING,
        OrderStatus.PROCESSING,
        OrderStatus.UPLOADING,
    ):
        ok = await runner.cancel_job(order_id)
        if ok:
            return {
                "status": "stopping",
                "order_id": order_id,
                "message": "SIGTERM sent to pipeline",
            }
        # 파이프라인은 이미 끝났고 annotation/업로드 중이면 PID가 없음 → 플래그만 남김
        runner.record_stop_request(order_id)
        return {
            "status": "cancelling",
            "order_id": order_id,
            "message": "Stop recorded; will cancel when the current step finishes (no pipeline process).",
        }

    raise HTTPException(
        status_code=400,
        detail=f"Order {order_id} is {job.status.value}; use Track for completed jobs",
    )


@app.post("/order/{order_id}/delete-run")
async def delete_order_run(order_id: str):
    """
    주문을 SQLite/메모리에서 제거하고, 해당 런의 analysis·output·log 디렉터리를 삭제합니다.
    FASTQ 디렉터리(`fastq_dir`)는 삭제하지 않습니다.
    실행 중(RUNNING 등)인 작업은 먼저 Stop 해야 합니다.
    메모리에 없고 SQLite에만 남은 꼬인 행은 자동으로 DB 정리(purge)를 시도합니다.
    """
    queue_manager = get_queue_manager()
    ok, message, detail = await queue_manager.delete_order_with_artifacts(order_id)
    if not ok:
        if "not found" in message.lower():
            ok_p, msg_p, det_p = await queue_manager.purge_order_db_only(order_id)
            if ok_p:
                det_p = dict(det_p)
                det_p["fallback"] = "db_only_not_in_memory"
                return {
                    "status": "ok",
                    "order_id": order_id,
                    "message": msg_p,
                    **det_p,
                }
        raise HTTPException(
            status_code=404 if "not found" in message.lower() else 400,
            detail=message,
        )
    return {"status": "ok", "order_id": order_id, "message": message, **detail}


@app.post("/order/{order_id}/purge-db")
async def purge_order_database_record(
    order_id: str,
    force: bool = Query(False),
):
    """
    SQLite orders 행 + 메모리/큐에서만 제거합니다. 디스크의 analysis/output/log 는 삭제하지 않습니다.
    job_json 파싱 실패 등으로 목록에 안 보이지만 DB에만 남은 경우, 스크립트로 같은 ID를 지우거나
    이 API를 호출하면 정리됩니다 (실행 중 작업은 Stop 후; stuck 시 force=true).
    """
    queue_manager = get_queue_manager()
    ok, message, detail = await queue_manager.purge_order_db_only(order_id, force=force)
    if not ok:
        raise HTTPException(status_code=400, detail=message)
    return {"status": "ok", "order_id": order_id, "message": message, **detail}


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

    job = _build_job_from_submit_request(service_code, request)

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


def _browse_fastq_directory_payload(path: str, service_code: Optional[str]) -> Dict[str, Any]:
    """
    Portal FASTQ 브라우저 공통 구현.
    루트: sgnipt → SGIPT_FASTQ_DIR, carrier_screening → CARRIER_SCREENING_FASTQ_DIR, 그 외 FASTQ_BASE_DIR.
    """
    root = _fastq_root_for_service(service_code)
    norm_prefix = _normalize_rel_path(path)

    if not norm_prefix and not os.path.isdir(root):
        env_hint = (
            "SGNIPT_FASTQ_DIR / CARRIER_SCREENING_FASTQ_DIR / FASTQ_BASE_DIR"
            if service_code
            else "FASTQ_BASE_DIR (or service-specific FASTQ dir)"
        )
        return {
            "root": root,
            "rel_path": "",
            "parent_rel": "",
            "root_exists": False,
            "service_code": service_code,
            "items": [],
            "hint": (
                "On this server, this service FASTQ root does not exist or is not a directory. "
                f"Create it or adjust {env_hint} in .env (resolved to: {root})."
            ),
        }

    target = _safe_join_under_fastq(path, service_code)
    if not os.path.exists(target):
        raise HTTPException(
            status_code=404,
            detail=f"Path not found under FASTQ root: {norm_prefix or '(root)'}",
        )
    if not os.path.isdir(target):
        raise HTTPException(status_code=400, detail="Not a directory")

    items: List[Dict[str, Any]] = []
    try:
        names = sorted(os.listdir(target))
    except OSError as e:
        raise HTTPException(status_code=403, detail=str(e)) from e
    for name in names:
        if name.startswith("."):
            continue
        full = os.path.join(target, name)
        rel_child = f"{norm_prefix}/{name}" if norm_prefix else name
        rel_child = rel_child.replace("\\", "/")
        if os.path.isdir(full):
            items.append({"name": name, "rel_path": rel_child, "kind": "dir"})
        elif os.path.isfile(full) and _is_fastq_filename(name):
            items.append({
                "name": name,
                "rel_path": rel_child,
                "abs_path": full,
                "kind": "file",
            })

    parent_rel = ""
    if norm_prefix:
        parent_rel = "/".join(norm_prefix.split("/")[:-1])

    return {
        "root": root,
        "rel_path": norm_prefix,
        "parent_rel": parent_rel,
        "root_exists": True,
        "service_code": service_code,
        "items": items,
    }


@app.get("/api/fastq/browse")
async def browse_fastq_short_path(
    path: str = Query(default="", description="서비스 FASTQ 루트 아래 상대 경로"),
    service_code: Optional[str] = Query(
        default=None,
        description="sgnipt | carrier_screening",
    ),
):
    """
    포털용 FASTQ 디렉터리 목록 (권장 경로).

    `/api/portal/browse/fastq` 와 동일. 일부 프록시/라우팅에서 `portal` 접두사가 막힐 때 사용.
    """
    return _browse_fastq_directory_payload(path, service_code)


@app.get("/api/portal/browse/fastq")
async def browse_fastq_directory(
    path: str = Query(default="", description="선택한 서비스 FASTQ 루트 아래 상대 경로"),
    service_code: Optional[str] = Query(
        default=None,
        description="sgnipt | carrier_screening (미지정 시 기본 fastq_base_dir)",
    ),
):
    """`/api/fastq/browse` 와 동일 (하위 호환)."""
    return _browse_fastq_directory_payload(path, service_code)


@app.patch("/order/{order_id}/fastq", response_model=None)
async def patch_order_fastq_paths(order_id: str, request: UpdateFastqPathsRequest):
    """
    QUEUED 주문의 로컬 R1/R2 경로만 수정 (실행 전).
    경로는 해당 주문 service_code의 FASTQ browse 루트 아래 파일이어야 합니다.
    """
    if request.fastq_r1_path is None and request.fastq_r2_path is None:
        raise HTTPException(status_code=400, detail="Provide fastq_r1_path and/or fastq_r2_path")

    new_r1: Optional[str] = None
    new_r2: Optional[str] = None
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    svc = job.service_code

    if request.fastq_r1_path is not None:
        s = request.fastq_r1_path.strip()
        new_r1 = "" if not s else _validate_optional_fastq_file(s, svc)
    if request.fastq_r2_path is not None:
        s = request.fastq_r2_path.strip()
        new_r2 = "" if not s else _validate_optional_fastq_file(s, svc)

    ok, err = await queue_manager.update_queued_job_fastq_paths(
        order_id,
        fastq_r1_path=new_r1,
        fastq_r2_path=new_r2,
    )
    if not ok:
        raise HTTPException(status_code=400, detail=err)

    job = queue_manager.get_job(order_id)
    return {
        "status": "ok",
        "order_id": order_id,
        "fastq_r1_path": job.fastq_r1_path if job else None,
        "fastq_r2_path": job.fastq_r2_path if job else None,
    }


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

    if job.service_code == "carrier_screening":
        from .services.carrier_screening.report import (
            carrier_report_template_kind,
            report_languages_from_order,
            _carrier_order_flat,
        )

        p_raw = job.params or {}
        kind = carrier_report_template_kind(p_raw)
        if kind is None:
            raise HTTPException(
                status_code=400,
                detail=(
                    "PDF report is only supported for standard carrier screening "
                    "or CouplesCarrier (Other test type)."
                ),
            )
        if report_languages_from_order(p_raw) is None:
            raise HTTPException(
                status_code=400,
                detail="Report language must be EN, CN, or KO for PDF generation.",
            )
        p = _carrier_order_flat(p_raw)
        if kind == "couples":
            p2 = (p.get("patient2_name") or "").strip()
            req_name = (
                (request.partner_info or {}).get("name") if request.partner_info else None
            )
            req_name = (str(req_name).strip() if req_name else "")
            if not p2 and not req_name:
                raise HTTPException(
                    status_code=400,
                    detail=(
                        "Couples report requires partner name "
                        "(order patient2_name or partner name in the form)."
                    ),
                )

    # 리포트 생성 (patient_info, partner_info; carrier_screening은 languages는 주문 Report Language로 결정)
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

    await queue_manager.mark_report_ready(
        order_id,
        message="Report ready for download",
    )

    store = queue_manager.store
    if store:
        review_langs = request.languages or []
        if job.service_code == "carrier_screening":
            from .services.carrier_screening.report import report_languages_from_order

            review_langs = report_languages_from_order(job.params or {}) or []
        review_blob = {
            "reviewer_info": request.reviewer_info,
            "patient_info": request.patient_info or {},
            "partner_info": request.partner_info,
            "languages": review_langs,
            "confirmed_variants": request.confirmed_variants,
            "saved_at": now_kst_iso(),
        }
        await asyncio.to_thread(store.set_review_json, order_id, review_blob)

    # 생성된 리포트 파일을 Platform에 업로드
    if job.service_code == "carrier_screening":
        from .services.carrier_screening.plugin import carrier_report_output_dir

        output_dir = carrier_report_output_dir(job)
    else:
        output_dir = job.output_dir or ""
    report_files = []

    if not output_dir:
        raise HTTPException(
            status_code=500,
            detail="Report output directory could not be resolved",
        )

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

    if store and output_dir:
        await asyncio.to_thread(
            ingest_report_json_from_disk, store, job, output_dir
        )

    # 플랫폼 업로드는 파일당 긴 타임아웃·순차 요청이라 포털이 멈춘 것처럼 보일 수 있음 → 백그라운드
    uploaded_count = 0
    platform_client = get_platform_client()
    if report_files and settings.platform_api_enabled:

        async def _upload_reports_background() -> None:
            try:
                upload_results = await platform_client.upload_all_outputs(
                    order_id, job.service_code, report_files
                )
                ok = sum(
                    1 for r in upload_results.values()
                    if r.status.value == "SUCCESS"
                )
                logger.info(
                    f"Background platform upload finished for {order_id}: "
                    f"{ok}/{len(report_files)} file(s)"
                )
            except Exception as e:
                logger.exception(
                    f"Background platform upload failed for {order_id}: {e}"
                )

        asyncio.create_task(_upload_reports_background())
        upload_msg = (
            f"Report generated: {len(report_files)} file(s). "
            "Platform upload is running in the background."
        )
    else:
        upload_msg = f"Report generated: {len(report_files)} file(s)."

    logger.info(
        f"Report generation complete for {order_id}: "
        f"{len(report_files)} files on disk (platform upload deferred={bool(report_files and settings.platform_api_enabled)})"
    )

    return ReportGenerateResponse(
        status="success",
        order_id=order_id,
        service_code=job.service_code,
        report_files=[os.path.basename(f.file_path) for f in report_files],
        uploaded_count=uploaded_count,
        message=upload_msg,
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
    for order_id, job in queue_manager._saved_jobs.items():
        all_jobs.append(job)
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
        # Portal "Completed" 필터가 분석 완료 후 리포트까지 포함하도록 REPORT_READY 동반 조회
        if status == OrderStatus.COMPLETED.value:
            all_jobs = [
                j for j in all_jobs
                if j.status.value in (OrderStatus.COMPLETED.value, OrderStatus.REPORT_READY.value)
            ]
        else:
            all_jobs = [j for j in all_jobs if j.status.value == status]

    # 최신순 정렬
    all_jobs.sort(key=lambda j: j.created_at or "", reverse=True)

    return {
        "orders": [
            {
                "order_id": j.order_id,
                "service_code": j.service_code,
                "sample_name": j.sample_name,
                "work_dir": j.work_dir,
                "status": j.status.value,
                "progress": j.progress,
                "message": j.message,
                "created_at": j.created_at,
                "started_at": j.started_at,
                "completed_at": j.completed_at,
                "updated_at": getattr(j, "updated_at", None),
                "priority": j.priority,
                "params": j.params or {},
                "fastq_r1_url": j.fastq_r1_url,
                "fastq_r2_url": j.fastq_r2_url,
                "fastq_r1_path": j.fastq_r1_path,
                "fastq_r2_path": j.fastq_r2_path,
                "callback_url": j.callback_url,
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
    디스크에 있으면 읽어 DB에도 스냅샷을 남기고, 없으면 SQLite에 저장된 스냅샷을 사용합니다.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    store = queue_manager.store
    result_json_path = None
    if job.output_dir:
        result_json_path = os.path.join(job.output_dir, "result.json")

    if result_json_path and os.path.isfile(result_json_path):
        with open(result_json_path, "r", encoding="utf-8") as f:
            result_data = json.load(f)
        if store and isinstance(result_data, dict):
            await asyncio.to_thread(store.set_result_json, order_id, result_data)
        if isinstance(result_data, dict):
            result_data = dict(result_data)
            result_data["order_params"] = job.params or {}
            return result_data
        return {"data": result_data, "order_params": job.params or {}}

    if store:
        cached = await asyncio.to_thread(store.get_result_json, order_id)
        if cached is not None:
            if isinstance(cached, dict):
                out = dict(cached)
                out["order_params"] = job.params or {}
                return out
            return {"data": cached, "order_params": job.params or {}}

    raise HTTPException(
        status_code=404,
        detail=(
            f"result.json not found for order: {order_id} "
            f"(no file under output_dir and no DB snapshot)."
        ),
    )


def _load_result_dict_for_gene_knowledge(order_id: str, job) -> Optional[Dict[str, Any]]:
    """Read result.json-like dict (variants list) without order_params merge."""
    result_json_path = None
    if job.output_dir:
        result_json_path = os.path.join(job.output_dir, "result.json")
    if result_json_path and os.path.isfile(result_json_path):
        with open(result_json_path, "r", encoding="utf-8") as f:
            return json.load(f)
    store = get_queue_manager().store
    if store:
        return store.get_result_json(order_id)
    return None


def _compute_order_gene_knowledge(
    order_id: str,
    job,
    enrich: bool,
    gene_filter: Optional[str] = None,
    force_refresh: bool = False,
    genes_csv: Optional[str] = None,
) -> Dict[str, Any]:
    """
    SQLite ``gene_data`` rows (Gemini cache) for each gene in the order's variant list.
    If ``enrich`` is True and ``GEMINI_API_KEY`` is set, call Gemini when narrative fields are
    missing (same extraction as report-time enrichment, including refresh when only disorder was cached).
    If ``gene_filter`` is set, only that symbol (must appear in the order's variants).
    If ``force_refresh`` is True, always re-run Gemini for the selected gene(s) (portal per-variant write-up).
    If ``genes_csv`` is not None, only genes in that comma-separated list (intersected with order variants);
    use ``""`` for an empty list (portal: no variants selected). If None, include all genes from the order.
    """
    from app.services.carrier_screening.gene_knowledge_db import (
        ensure_gene_knowledge_full_text,
        init_gene_knowledge_database,
        load_variant_knowledge_for_keys,
        make_variant_key,
        read_gene_knowledge_full_row,
        refresh_gene_knowledge_from_gemini,
    )

    db_path = (settings.gene_knowledge_db or "").strip()
    if not db_path:
        msg = "gene_knowledge_db is not configured on the daemon"
        return {
            "gene_knowledge_db_configured": False,
            "genes": {},
            "variants": {},
            "message": msg,
            "error": msg,
        }

    init_gene_knowledge_database(db_path)
    result_data = _load_result_dict_for_gene_knowledge(order_id, job)
    if not isinstance(result_data, dict):
        msg = "result.json not available for this order"
        return {
            "gene_knowledge_db_configured": True,
            "genes": {},
            "variants": {},
            "message": msg,
            "error": msg,
        }

    variants = result_data.get("variants") or []
    order_genes = {
        (v.get("gene") or "").strip().upper() for v in variants if (v.get("gene") or "").strip()
    }
    genes = sorted(order_genes)
    gf = (gene_filter or "").strip().upper() or None
    if gf:
        if gf not in order_genes:
            return {
                "gene_knowledge_db_configured": True,
                "genes": {},
                "variants": {},
                "error": f"Gene {gf} is not among this order's variants",
                "enrich_requested": enrich,
                "force_refresh": force_refresh,
                "gemini_available": bool((settings.gemini_api_key or "").strip()),
            }
        genes = [gf]
    elif genes_csv is not None:
        raw = (genes_csv or "").strip()
        if not raw:
            genes = []
        else:
            requested = {x.strip().upper() for x in raw.split(",") if x.strip()}
            genes = sorted(order_genes & requested)

    gemini_key = (settings.gemini_api_key or "").strip()
    model = getattr(settings, "gene_knowledge_gemini_model", "gemini-2.5-flash")
    allow_gemini = bool(enrich and gemini_key)

    gene_set = set(genes)
    variant_keys_set: set = set()
    for v in variants:
        g = (v.get("gene") or "").strip().upper()
        if not g or g not in gene_set:
            continue
        variant_keys_set.add(
            make_variant_key(g, str(v.get("hgvsc") or ""), str(v.get("hgvsp") or ""))
        )
    variants_out = load_variant_knowledge_for_keys(db_path, list(variant_keys_set))

    if force_refresh and not gemini_key:
        return {
            "gene_knowledge_db_configured": True,
            "genes": {},
            "variants": variants_out,
            "error": "GEMINI_API_KEY not configured on the daemon",
            "enrich_requested": enrich,
            "force_refresh": True,
            "gemini_available": False,
        }

    out: Dict[str, Any] = {}
    gemini_fetch_error: Optional[str] = None
    for gene in genes:
        if force_refresh and gemini_key:
            row, gerr = refresh_gene_knowledge_from_gemini(gene, db_path, gemini_key, model=model)
            if gerr:
                gemini_fetch_error = gerr
        elif allow_gemini:
            row = ensure_gene_knowledge_full_text(
                gene, db_path, gemini_key, model=model, allow_gemini=True
            )
        else:
            row = read_gene_knowledge_full_row(gene, db_path)
        out[gene] = row if row else {}

    ret: Dict[str, Any] = {
        "gene_knowledge_db_configured": True,
        "genes": out,
        "variants": variants_out,
        "enrich_requested": enrich,
        "force_refresh": force_refresh,
        "gemini_available": bool(gemini_key),
    }
    if force_refresh and gemini_fetch_error:
        ret["gemini_fetch_error"] = gemini_fetch_error
    return ret


def _put_order_gene_knowledge(
    order_id: str,
    job,
    body: GeneKnowledgeSaveRequest,
    genes_csv: Optional[str] = None,
) -> Dict[str, Any]:
    """Write one ``gene_data`` row; gene must be in this order's variant list (and in ``genes_csv`` if given)."""
    from app.services.carrier_screening.gene_knowledge_db import (
        init_gene_knowledge_database,
        read_gene_knowledge_full_row,
        upsert_gene_data,
    )

    db_path = (settings.gene_knowledge_db or "").strip()
    if not db_path:
        raise HTTPException(
            status_code=503,
            detail="gene_knowledge_db is not configured on the daemon",
        )

    gene = (body.gene or "").strip().upper()
    if not gene:
        raise HTTPException(status_code=400, detail="gene is required")

    result_data = _load_result_dict_for_gene_knowledge(order_id, job)
    if not isinstance(result_data, dict):
        raise HTTPException(
            status_code=400,
            detail="result.json not available for this order",
        )

    variants = result_data.get("variants") or []
    order_genes = {
        (v.get("gene") or "").strip().upper()
        for v in variants
        if (v.get("gene") or "").strip()
    }
    if gene not in order_genes:
        raise HTTPException(
            status_code=400,
            detail=f"Gene {gene} is not among this order's variants",
        )

    if genes_csv is not None:
        raw = (genes_csv or "").strip()
        allowed = {x.strip().upper() for x in raw.split(",") if x.strip()} if raw else set()
        if gene not in allowed:
            raise HTTPException(
                status_code=400,
                detail="Gene is not in the current selection (genes query does not include this symbol)",
            )

    init_gene_knowledge_database(db_path)
    upsert_gene_data(
        db_path,
        {
            "gene_symbol": gene,
            "function_summary": body.function_summary or "",
            "disease_association": body.disease_association or "",
            "disorder": body.disorder or "",
            "omim_number": body.omim_number or "",
            "inheritance": body.inheritance or "",
        },
    )
    row = read_gene_knowledge_full_row(gene, db_path)
    return {"ok": True, "gene": gene, "row": row or {}}


def _put_order_variant_knowledge(
    order_id: str,
    job,
    body: VariantKnowledgeSaveRequest,
    genes_csv: Optional[str] = None,
) -> Dict[str, Any]:
    """Upsert ``variant_knowledge`` for one variant (stable key across orders)."""
    from app.services.carrier_screening.gene_knowledge_db import (
        init_gene_knowledge_database,
        make_variant_key,
        read_variant_knowledge_row,
        upsert_variant_knowledge,
    )

    db_path = (settings.gene_knowledge_db or "").strip()
    if not db_path:
        raise HTTPException(
            status_code=503,
            detail="gene_knowledge_db is not configured on the daemon",
        )

    vk = (body.variant_key or "").strip()
    if not vk:
        raise HTTPException(status_code=400, detail="variant_key is required")

    result_data = _load_result_dict_for_gene_knowledge(order_id, job)
    if not isinstance(result_data, dict):
        raise HTTPException(
            status_code=400,
            detail="result.json not available for this order",
        )

    variants = result_data.get("variants") or []
    matched: Optional[Dict[str, Any]] = None
    for v in variants:
        g = (v.get("gene") or "").strip().upper()
        if not g:
            continue
        if make_variant_key(g, str(v.get("hgvsc") or ""), str(v.get("hgvsp") or "")) == vk:
            matched = v
            break
    if not matched:
        raise HTTPException(
            status_code=400,
            detail="variant_key does not match any variant in this order",
        )

    gene = (matched.get("gene") or "").strip().upper()
    if genes_csv is not None:
        raw = (genes_csv or "").strip()
        allowed = {x.strip().upper() for x in raw.split(",") if x.strip()} if raw else set()
        if gene not in allowed:
            raise HTTPException(
                status_code=400,
                detail="Gene is not in the current selection (genes query does not include this symbol)",
            )

    init_gene_knowledge_database(db_path)
    upsert_variant_knowledge(
        db_path,
        {
            "variant_key": vk,
            "gene_symbol": gene,
            "hgvsc": str(matched.get("hgvsc") or ""),
            "hgvsp": str(matched.get("hgvsp") or ""),
            "variant_notes": body.variant_notes or "",
        },
    )
    row = read_variant_knowledge_row(vk, db_path)
    return {"ok": True, "variant_key": vk, "row": row or {}}


@app.put("/order/{order_id}/variant-knowledge")
async def put_order_variant_knowledge(
    order_id: str,
    body: VariantKnowledgeSaveRequest,
    genes: Optional[str] = Query(
        None,
        description="Comma-separated gene symbols allowed (portal: selected variants)",
    ),
):
    """Save per-variant notes to ``variant_knowledge`` (gene-level text remains in ``gene_data``)."""
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    return await asyncio.to_thread(_put_order_variant_knowledge, order_id, job, body, genes)


@app.put("/order/{order_id}/gene-knowledge")
async def put_order_gene_knowledge(
    order_id: str,
    body: GeneKnowledgeSaveRequest,
    genes: Optional[str] = Query(
        None,
        description="Comma-separated gene symbols allowed (portal: selected variants); if set, body.gene must be listed",
    ),
):
    """Save portal edits to the shared ``gene_knowledge`` SQLite cache for one gene."""
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    return await asyncio.to_thread(_put_order_gene_knowledge, order_id, job, body, genes)


@app.get("/order/{order_id}/gene-knowledge")
async def get_order_gene_knowledge(
    order_id: str,
    enrich: bool = Query(False, description="If true, call Gemini for genes missing from SQLite cache"),
    gene: Optional[str] = Query(
        None,
        description="If set, only this gene symbol (must appear in order variants)",
    ),
    force: bool = Query(
        False,
        description="If true, always re-run Gemini for the selected gene(s) (new write-up)",
    ),
    genes: Optional[str] = Query(
        None,
        description=(
            "Comma-separated gene symbols to include (must each appear in order variants). "
            "If omitted, all order genes are included. Use an empty value (genes=) for none."
        ),
    ),
):
    """
    Gene-level text from the ``gene_knowledge`` SQLite cache (same source as report-time enrichment:
    ``function_summary``, ``disease_association``, disorder, OMIM, inheritance).

    Use this to populate the portal Gene database tab and richer Review Case report narratives.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    return await asyncio.to_thread(
        _compute_order_gene_knowledge, order_id, job, enrich, gene, force, genes
    )


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


def _safe_order_file_path(root: str, rel: str) -> Optional[str]:
    """rel 은 주문 출력 루트 아래 상대 경로 (예 qc/plot.png). 경로 탈출 방지."""
    if not rel or not root or ".." in rel.replace("\\", "/"):
        return None
    rel = rel.strip().lstrip("/").replace("\\", "/")
    try:
        root_abs = os.path.realpath(root)
        full = os.path.realpath(os.path.join(root, *rel.split("/")))
    except OSError:
        return None
    if full != root_abs and not full.startswith(root_abs + os.sep):
        return None
    return full if os.path.isfile(full) else None


def _resolve_order_artifact_path(job: Job, filename: str) -> Optional[str]:
    """
    artifact output_dir 우선, carrier_screening 은 layout/script output 에 동일 상대 경로 시도.
    """
    rel = (filename or "").strip().lstrip("/")
    if not rel or ".." in rel.replace("\\", "/"):
        return None
    roots: List[str] = []
    if job.output_dir:
        roots.append(job.output_dir)
    if (job.service_code or "").strip() == "carrier_screening":
        from .services.carrier_screening.layout_norm import carrier_sequencing_folder

        wk = str(job.work_dir).strip() or "00"
        smp = str(job.sample_name).strip()
        seq = carrier_sequencing_folder(job)
        layout_base = (settings.carrier_screening_layout_base or "").strip()
        if layout_base and smp:
            roots.append(os.path.join(layout_base, "output", wk, smp))
        if seq and seq != smp and layout_base:
            roots.append(os.path.join(layout_base, "output", wk, seq))
        sd = (settings.carrier_screening_script_data_dir or "").strip()
        if sd and smp:
            roots.append(os.path.join(sd, "output", wk, smp))
        if sd and seq and seq != smp:
            roots.append(os.path.join(sd, "output", wk, seq))

    seen: set = set()
    for root in roots:
        if not root or not os.path.isdir(root):
            continue
        try:
            key = os.path.realpath(root)
        except OSError:
            continue
        if key in seen:
            continue
        seen.add(key)
        hit = _safe_order_file_path(root, rel)
        if hit:
            return hit
    return None


@app.get("/order/{order_id}/file/{filename:path}")
async def download_order_file(order_id: str, filename: str):
    """
    주문의 특정 출력 파일을 다운로드합니다.
    filename 에 슬래시 포함 가능 (예 qc/plot.png).
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    file_path = _resolve_order_artifact_path(job, filename)
    if not file_path:
        raise HTTPException(status_code=404, detail=f"File not found: {filename}")

    base_name = os.path.basename(filename) or "download"
    return FileResponse(
        path=file_path,
        filename=base_name,
        media_type=_guess_content_type(filename),
    )


@app.get("/order/{order_id}/pipeline-log")
async def get_order_pipeline_log(
    order_id: str,
    max_bytes: int = Query(
        default=524_288,
        ge=4096,
        le=4_194_304,
        description="큰 로그 파일은 끝에서부터 이 크기(바이트)만 반환",
    ),
):
    """
    `log_dir/pipeline.log` 내용을 plain text로 반환합니다 (Portal 주문 상세 터미널).
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    log_dir = job.log_dir
    if not log_dir:
        raise HTTPException(status_code=404, detail="log_dir not set for this order")
    log_dir_real = os.path.realpath(log_dir)
    candidate = os.path.join(log_dir, "pipeline.log")
    file_path = os.path.realpath(candidate)
    if not file_path.startswith(log_dir_real + os.sep):
        raise HTTPException(status_code=400, detail="invalid log path")
    if not os.path.isfile(file_path):
        raise HTTPException(
            status_code=404,
            detail="pipeline.log not found (pipeline may not have started yet)",
        )
    try:
        size = os.path.getsize(file_path)
        with open(file_path, "rb") as f:
            if size > max_bytes:
                f.seek(size - max_bytes)
                f.readline()
                raw = f.read()
            else:
                raw = f.read()
        text = raw.decode("utf-8", errors="replace")
    except OSError as e:
        logger.warning("pipeline-log read failed for %s: %s", order_id, e)
        raise HTTPException(status_code=500, detail=f"cannot read pipeline.log: {e}")
    return PlainTextResponse(content=text, media_type="text/plain; charset=utf-8")


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
        "timestamp": now_kst_iso()
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
    job.created_at = now_kst_iso()
    job.progress = 100 if status_enum == OrderStatus.COMPLETED else 0

    if status_enum == OrderStatus.COMPLETED:
        job.completed_at = now_kst_iso()
        job.started_at = now_kst_iso()
        job.message = "Mock job - analysis complete"
        queue_manager._completed_jobs[order_id] = job
    elif status_enum == OrderStatus.RUNNING:
        job.started_at = now_kst_iso()
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


# ═══════════════════════════════════════════════════════════════
# Literature API
# ─────────────────────────────────────────────────────────────
# 모든 포털에서 호출 가능한 안정적인 REST API.
# Portal(index.html)은 이 API를 소비하는 얇은 클라이언트.
# NIPT Client Portal 이식 시에도 이 엔드포인트를 그대로 사용.
# ═══════════════════════════════════════════════════════════════

@app.get("/api/literature/search")
async def literature_search(
    gene: str = Query(..., description="유전자 심볼 (예: BRCA2)"),
    hgvsc: str = Query("", description="HGVS.c 표기 (예: c.5266dupC)"),
    hgvsp: str = Query("", description="HGVS.p 표기 (예: p.Gln1756fs)"),
    effect: str = Query("", description="변이 효과 (예: frameshift_variant)"),
    force_refresh: bool = Query(False, description="캐시 무시하고 재검색"),
    max_results: int = Query(10, ge=1, le=30),
):
    """
    변이에 대한 PubMed 문헌을 검색합니다. 결과는 영구 캐시됩니다.

    어떤 포털에서도 호출 가능한 안정적인 API 계약:
      GET /api/literature/search?gene=BRCA2&hgvsc=c.5266dupC
    """
    if not settings.literature_enabled:
        raise HTTPException(status_code=503, detail="Literature search is disabled (LITERATURE_ENABLED=false)")

    from .services.carrier_screening.literature import search_variant_literature
    result = await search_variant_literature(
        gene=gene,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        effect=effect,
        max_results=max_results,
        ncbi_email=settings.ncbi_email,
        ncbi_api_key=settings.ncbi_api_key,
        ncbi_tool=settings.ncbi_tool,
        force_refresh=force_refresh,
    )
    return result


@app.get("/api/literature/articles")
def literature_list_articles(
    cursor: int = Query(0, ge=0),
    count: int = Query(50, ge=1, le=200),
    search: str = Query("", description="제목/초록/저자 검색"),
    sort_by: str = Query("cached_at", description="정렬 기준: cached_at | pub_date"),
):
    """
    캐시된 PubMed 논문 목록을 반환합니다.
    포털의 '문헌 관리' 페이지에서 사용합니다.
    """
    if not settings.literature_enabled:
        raise HTTPException(status_code=503, detail="Literature search is disabled")

    from .services.carrier_screening.literature import list_cached_articles
    return list_cached_articles(cursor=cursor, count=count, search=search, sort_by=sort_by)


@app.get("/api/literature/articles/stats")
def literature_stats():
    """캐시 통계 (논문 수, 유전자 수 등)."""
    if not settings.literature_enabled:
        return {"enabled": False}

    from .services.carrier_screening.literature import get_literature_stats
    stats = get_literature_stats()
    stats["enabled"] = True
    return stats


@app.get("/api/literature/articles/{pmid}")
def literature_get_article(pmid: str):
    """PMID로 캐시된 논문 상세 조회."""
    from .services.carrier_screening.literature import get_cached_article
    art = get_cached_article(pmid)
    if not art:
        raise HTTPException(status_code=404, detail=f"Article {pmid} not found in cache")
    return art


@app.delete("/api/literature/articles/{pmid}")
def literature_delete_article(pmid: str):
    """캐시에서 논문 1건 삭제."""
    from .services.carrier_screening.literature import delete_cached_article
    deleted = delete_cached_article(pmid)
    return {"deleted": deleted, "pmid": pmid}


@app.delete("/api/literature/cache")
def literature_clear_cache():
    """전체 캐시 삭제 (관리자 기능)."""
    from .services.carrier_screening.literature import clear_literature_cache
    result = clear_literature_cache()
    return result


@app.get("/api/literature/searches")
def literature_list_searches(
    gene: str = Query("", description="유전자 필터"),
    cursor: int = Query(0, ge=0),
    count: int = Query(50, ge=1, le=200),
):
    """변이별 검색 캐시 목록 (어떤 변이가 검색되었는지 확인)."""
    from .services.carrier_screening.literature import list_cached_searches
    return list_cached_searches(gene=gene, cursor=cursor, count=count)
