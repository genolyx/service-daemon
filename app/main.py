"""
Service Daemon - FastAPI Application

여러 유전체 분석 서비스를 통합 관리하는 범용 데몬입니다.
기존 nipt-daemon의 구조를 일반화하여 플러그인 기반으로 동작합니다.
"""

import os
import re
import json
import asyncio
import glob
import hmac
import logging
from urllib.parse import unquote, urlparse
from contextlib import asynccontextmanager
from typing import Dict, Any, List, Optional, Tuple

from fastapi import FastAPI, HTTPException, Query, Body, Request
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse, FileResponse, PlainTextResponse, RedirectResponse
from fastapi.middleware.cors import CORSMiddleware
from uvicorn.middleware.proxy_headers import ProxyHeadersMiddleware

from .config import settings
from .datetime_kst import now_kst_iso, now_kst_date_compact
from .logging_config import setup_logging, setup_middleware
from .models import (
    OrderSubmitRequest, OrderSubmitResponse, OrderSaveResponse, OrderStatusResponse,
    OrderUpdateRequest, OrderUpdateResponse, StartOrderRequest,
    OrderStatus, Job, QueueSummary, OutputFile,
    ReportGenerateRequest, ReportGenerateResponse,
    GeneKnowledgeSaveRequest,
    VariantKnowledgeSaveRequest,
    UpdateFastqPathsRequest,
    DarkGenesReviewRequest,
    PgxReviewRequest,
    WesPanelCustomSave,
)
from .queue_manager import get_queue_manager
from .order_store import ingest_report_json_from_disk
from .annotation_resources import annotation_resource_report
from .runner import get_runner
from .platform_client import get_platform_client
from .services import load_plugins, get_plugin, list_service_codes, get_all_plugins
from .services.carrier_screening.prior_reuse import prior_reuse_artifact_roots

logger = logging.getLogger(__name__)

_FASTQ_NAME_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# Portal GET /order/{id}/result must never be served from browser/CDN cache (stale variants / dark_genes).
_RESULT_JSON_CACHE_HEADERS = {
    "Cache-Control": "no-store, no-cache, must-revalidate",
    "Pragma": "no-cache",
}

# Same Nextflow + carrier_screening review/report stack (optional WES panel for whole_exome).
_CARRIER_LIKE = frozenset({"carrier_screening", "whole_exome", "health_screening"})


def _order_result_review_debug(
    job: Job,
    result_json_path: Optional[str],
    source: str,
    result_payload: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Operator-facing paths for diagnosing shared result.json between orders.
    Returned only when GET /order/{id}/result?debug=1.
    """
    out: Dict[str, Any] = {
        "order_id": job.order_id,
        "service_code": job.service_code,
        "work_dir": getattr(job, "work_dir", None),
        "sample_folder": getattr(job, "sample_name", None),
        "job_output_dir": getattr(job, "output_dir", None),
        "job_analysis_dir": getattr(job, "analysis_dir", None),
        "resolved_result_json_path": result_json_path,
        "resolved_path_exists": bool(result_json_path and os.path.isfile(result_json_path)),
        "source": source,
    }
    if job.service_code in _CARRIER_LIKE:
        from app.services.carrier_screening.plugin import carrier_report_output_dir
        from app.services.wes_panels import (
            get_panel_by_id,
            interpretation_gene_set_for_job,
            resolve_panel_interpretation_genes,
            should_apply_interpretation_post_filter,
        )

        try:
            cro = carrier_report_output_dir(job)
            out["carrier_report_output_dir"] = cro
            canon = os.path.join(cro, "result.json")
            out["canonical_result_json_path"] = canon
            out["canonical_path_exists"] = os.path.isfile(canon)
        except Exception as e:
            out["carrier_report_output_dir_error"] = str(e)

        p = job.params or {}
        wpid = p.get("wes_panel_id")
        if not wpid and isinstance(p.get("carrier"), dict):
            wpid = p["carrier"].get("wes_panel_id")
        wpid_s = (str(wpid).strip() if wpid else None)
        out["wes_panel_id"] = wpid_s
        if wpid_s:
            panel = get_panel_by_id(wpid_s)
            out["panel_found_in_catalog"] = panel is not None
            if panel:
                out["panel_catalog_interpretation_gene_count"] = len(
                    resolve_panel_interpretation_genes(panel)
                )
        pig = p.get("panel_interpretation_genes")
        out["job_panel_interpretation_genes_count"] = (
            len(pig) if isinstance(pig, list) else None
        )
        out["order_interpretation_gene_count"] = len(interpretation_gene_set_for_job(job))
        out["interpretation_post_filter_enabled"] = should_apply_interpretation_post_filter(job)

        pr = (job.params or {}).get("_prior_reuse")
        if pr:
            out["prior_pipeline_reuse"] = True
            out["prior_reuse_order_id"] = (job.params or {}).get("_prior_reuse_order_id")

    if isinstance(result_payload, dict):
        fs = result_payload.get("filter_summary") or {}
        if isinstance(fs, dict) and fs:
            out["result_filter_summary"] = {
                "wes_panel_id": fs.get("wes_panel_id"),
                "wes_panel_label": fs.get("wes_panel_label"),
                "interpretation_post_filter_genes": fs.get("interpretation_post_filter_genes"),
                "interpretation_post_filter_applied": fs.get("interpretation_post_filter_applied"),
            }
    return out


def _fastq_root_real() -> str:
    return os.path.realpath(settings.fastq_base_dir)


def _fastq_root_for_service(service_code: Optional[str]) -> str:
    """포털 browse / FASTQ 경로 검증용 루트 (서비스별)."""
    sc = (service_code or "").strip().lower().replace("-", "_")
    if sc == "sgnipt":
        return os.path.realpath(settings.sgnipt_fastq_root)
    if sc in ("carrier_screening", "whole_exome", "health_screening"):
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


def _is_csv_filename(name: str) -> bool:
    return name.lower().endswith(".csv")


def _bam_csv_root_for_service(service_code: Optional[str]) -> str:
    """BAM samplesheet CSV 브라우져용 루트 디렉터리 (서비스별)."""
    sc = (service_code or "").strip().lower().replace("-", "_")
    if sc in ("carrier_screening", "whole_exome", "health_screening"):
        return os.path.realpath(os.path.join(settings.carrier_screening_work_root, "data"))
    # sgnipt (default): SGNIPT_DATA_DIR 우선, 없으면 sgnipt_job_root/data
    data_dir = (settings.sgnipt_data_dir or "").strip()
    if data_dir:
        return os.path.realpath(data_dir)
    return os.path.realpath(os.path.join(settings.sgnipt_job_root, "data"))


def _browse_bam_csv_directory_payload(
    path: str,
    service_code: Optional[str],
    abs_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    BAM samplesheet CSV 브라우져 공통 구현.

    abs_path 가 주어지면 임의 절대 경로를 직접 탐색 (루트 제한 없음).
    path 는 서비스별 데이터 루트 기준 상대 경로.
    """
    default_root = _bam_csv_root_for_service(service_code)

    # abs_path 직접 이동 모드
    if abs_path:
        target = os.path.realpath(abs_path.strip())
        if not os.path.exists(target):
            raise HTTPException(status_code=404, detail=f"Path not found: {abs_path}")
        if not os.path.isdir(target):
            raise HTTPException(status_code=400, detail=f"Not a directory: {abs_path}")
        root = target  # 이 디렉토리를 임시 루트로 사용
        norm_prefix = ""
        parent_abs = os.path.dirname(target)
    else:
        root = default_root
        norm_prefix = _normalize_rel_path(path)

        if not norm_prefix and not os.path.isdir(root):
            return {
                "root": root,
                "rel_path": "",
                "parent_rel": "",
                "root_exists": False,
                "service_code": service_code,
                "items": [],
                "hint": (
                    f"BAM data root does not exist: {root}. "
                    "Create the directory or check SGNIPT_DATA_DIR / SGNIPT_WORK_ROOT in .env."
                ),
            }

        target = os.path.realpath(os.path.join(root, norm_prefix) if norm_prefix else root)
        root_sep = root if root.endswith(os.sep) else root + os.sep
        if target != root and not target.startswith(root_sep):
            raise HTTPException(status_code=400, detail="Path escapes BAM data root")

        if not os.path.exists(target):
            raise HTTPException(
                status_code=404,
                detail=f"Path not found under BAM data root: {norm_prefix or '(root)'}",
            )
        if not os.path.isdir(target):
            raise HTTPException(status_code=400, detail="Not a directory")

        parent_abs = os.path.dirname(target) if norm_prefix else None

    items: List[Dict[str, Any]] = []
    try:
        names = sorted(os.listdir(target))
    except OSError as e:
        raise HTTPException(status_code=403, detail=str(e)) from e

    for name in names:
        if name.startswith("."):
            continue
        full = os.path.join(target, name)
        if abs_path:
            rel_child = name
            child_abs = full
        else:
            rel_child = f"{norm_prefix}/{name}" if norm_prefix else name
            rel_child = rel_child.replace("\\", "/")
            child_abs = full
        if os.path.isdir(full):
            entry: Dict[str, Any] = {"name": name, "kind": "dir"}
            if abs_path:
                entry["abs_path"] = child_abs
            else:
                entry["rel_path"] = rel_child
            items.append(entry)
        elif os.path.isfile(full) and _is_csv_filename(name):
            entry = {"name": name, "abs_path": child_abs, "kind": "file"}
            if not abs_path:
                entry["rel_path"] = rel_child
            items.append(entry)

    if abs_path:
        return {
            "root": root,
            "current_abs": target,
            "parent_abs": parent_abs,
            "rel_path": "",
            "parent_rel": "",
            "root_exists": True,
            "service_code": service_code,
            "items": items,
        }

    parent_rel = ""
    if norm_prefix:
        parent_rel = "/".join(norm_prefix.split("/")[:-1])

    return {
        "root": root,
        "current_abs": target,
        "rel_path": norm_prefix,
        "parent_rel": parent_rel,
        "root_exists": True,
        "service_code": service_code,
        "items": items,
    }


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

_OPENAPI_DESCRIPTION = """
## Portal (UI)

로컬·개발 시 브라우저에서 **`/portal/`** 로 접속합니다. (예: `http://localhost:8003/portal/`)

- 주문·큐·Review 등은 이 UI가 동일 오리진의 **`/order/...`** API를 호출합니다.
- 리버스 프록시가 **`/api/order/...`** 만 넘기는 경우에도, 미들웨어가 내부적으로 `/order/...` 로 변환합니다.

## API 문서 (이 페이지)

- **Swagger UI**: `/docs` (현재 페이지)
- **ReDoc**: `/redoc`
- **OpenAPI JSON**: `/openapi.json`

## 선택적 인증

환경 변수 **`API_KEY`** 가 설정되어 있으면, 대부분의 API는 **`Authorization: Bearer <API_KEY>`** 또는 **`X-API-Key: <API_KEY>`** 가 필요합니다.

예외(키 없이 접근 가능): `/health`, `/portal/*`, `/docs`, `/openapi.json`, `/redoc`, 일부 `/api/portal/*` 등 — 상세는 `api_access_key_guard` 참고.

## 개발 문서

레포지토리 **`docs/API_DEVELOPMENT.md`** 에 Proxy·Portal·curl 예시를 정리해 두었습니다.
""".strip()

app = FastAPI(
    title="Service Daemon",
    description=_OPENAPI_DESCRIPTION,
    version="2.0.0",
    lifespan=lifespan,
)

# CORS 미들웨어 (테스트 Portal에서 로컬 API 호출 허용)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
# Behind nginx / SSH tunnels: honor X-Forwarded-* so request.url reflects the client-facing host/scheme.
app.add_middleware(ProxyHeadersMiddleware)

# 미들웨어 설정
setup_middleware(app)


def _strip_api_prefix_for_order_routes(path: str) -> str:
    """
    Reverse proxies sometimes only forward ``/api/*`` to this app. Expose the same
    order/queue/health routes under ``/api/order/...`` etc. by rewriting to the real paths.
    Do not strip ``/api/literature``, ``/api/fastq``, ``/api/portal`` — those are native.
    """
    if not path.startswith("/api/"):
        return path
    if path.startswith("/api/literature") or path.startswith("/api/fastq") or path.startswith("/api/portal"):
        return path
    if path == "/api/health":
        return "/health"
    for prefix in ("/api/order", "/api/orders", "/api/queue", "/api/services"):
        if path == prefix or path.startswith(prefix + "/"):
            return path[4:]
    return path


@app.middleware("http")
async def strip_api_prefix_middleware(request: Request, call_next):
    """Run before api_access_key_guard so paths match exemptions after rewrite."""
    path = request.url.path
    new_path = _strip_api_prefix_for_order_routes(path)
    if new_path != path:
        request.scope["path"] = new_path
        request.scope["raw_path"] = new_path.encode("utf-8")
    return await call_next(request)


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
    if (
        path == "/health"
        or path.startswith("/portal")
        or path.startswith("/api/portal")
        or path.startswith("/api/wes-panels")
        or path.startswith("/static")
        # OpenAPI / Swagger — 개발 시 API 탐색용 (API_KEY 있어도 문서 접근 가능)
        or path.startswith("/docs")
        or path == "/openapi.json"
        or path.startswith("/redoc")
    ):
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
            "whole_exome": _fastq_root_for_service("whole_exome"),
            "health_screening": _fastq_root_for_service("health_screening"),
            "default": fq_default,
        },
    }
    if annotation:
        body["annotation_resources"] = annotation_resource_report()
    return body


# ─── Order Endpoints ──────────────────────────────────────

def _pipeline_folder_label(order_id: str, sample_name: Optional[str]) -> str:
    """work_dir 아래 디렉터리 명: sample_name 비어 있으면 order_id."""
    sn = (sample_name or "").strip()
    oid = (order_id or "").strip()
    return sn or oid


def _build_job_from_submit_request(service_code: str, request: OrderSubmitRequest) -> Job:
    """Submit/Save 공통 Job 생성."""
    work_dir = request.work_dir or now_kst_date_compact()
    folder = _pipeline_folder_label(request.order_id, request.sample_name)
    if service_code == "sgnipt":
        root = settings.sgnipt_job_root
        oid = (request.order_id or "").strip()
        return Job(
            order_id=request.order_id,
            service_code=service_code,
            sample_name=folder,
            work_dir=work_dir,
            fastq_r1_url=request.fastq_r1_url,
            fastq_r2_url=request.fastq_r2_url,
            fastq_r1_path=request.fastq_r1_path,
            fastq_r2_path=request.fastq_r2_path,
            params=request.params or {},
            fastq_dir=os.path.join(root, "fastq", work_dir, oid),
            analysis_dir=os.path.join(root, "analysis", work_dir, oid),
            output_dir=os.path.join(root, "output", work_dir, oid),
            log_dir=os.path.join(root, "log", work_dir, oid),
        )
    if service_code in ("carrier_screening", "whole_exome", "health_screening"):
        work_root = settings.carrier_screening_work_root
        return Job(
            order_id=request.order_id,
            service_code=service_code,
            sample_name=folder,
            work_dir=work_dir,
            fastq_r1_url=request.fastq_r1_url,
            fastq_r2_url=request.fastq_r2_url,
            fastq_r1_path=request.fastq_r1_path,
            fastq_r2_path=request.fastq_r2_path,
            params=request.params or {},
            fastq_dir=os.path.join(work_root, "fastq", work_dir, folder),
            analysis_dir=os.path.join(work_root, "analysis", work_dir, folder),
            output_dir=os.path.join(work_root, "output", work_dir, folder),
            log_dir=os.path.join(work_root, "log", work_dir, folder),
        )
    base = settings.base_dir
    return Job(
        order_id=request.order_id,
        service_code=service_code,
        sample_name=folder,
        work_dir=work_dir,
        fastq_r1_url=request.fastq_r1_url,
        fastq_r2_url=request.fastq_r2_url,
        fastq_r1_path=request.fastq_r1_path,
        fastq_r2_path=request.fastq_r2_path,
            params=request.params or {},
            fastq_dir=os.path.join(base, "fastq", work_dir, folder),
        analysis_dir=os.path.join(base, "analysis", work_dir, folder),
        output_dir=os.path.join(base, "output", work_dir, folder),
        log_dir=os.path.join(base, "log", work_dir, folder),
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
    # Save = 초안 저장: wes_panel_id 등 Submit 전에 채울 수 있는 항목은 strict=False 로 건너뜀
    is_valid, error_msg = plugin.validate_params(request.params or {}, strict=False)
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
async def start_saved_order(
    order_id: str,
    request: Optional[StartOrderRequest] = Body(default=None),
):
    """
    SAVED 주문을 큐에 넣거나, FAILED/CANCELLED/COMPLETED 주문을 재시도로 다시 큐에 넣습니다.

    - body 없음 또는 ``{"fresh": false}`` (기본): Nextflow -resume 등 캐시 재활용 (Force Run).
    - ``{"fresh": true}``: 파이프라인 캐시를 삭제하고 처음부터 재실행 (Force Run Fresh / ``--fresh``).
    """
    fresh = False
    use_ssd = False
    scratch_dir = None
    if request is not None:
        fresh = request.fresh
        use_ssd = request.use_ssd
        scratch_dir = request.scratch_dir
    queue_manager = get_queue_manager()
    job_preview = queue_manager.get_job(order_id)
    if job_preview and job_preview.service_code in _CARRIER_LIKE:
        pl = get_plugin(job_preview.service_code)
        if pl:
            ok, err = pl.validate_params(job_preview.params or {})
            if not ok:
                raise HTTPException(status_code=400, detail=f"Invalid params: {err}")
    try:
        job, queue_position = await queue_manager.start_saved_job(
            order_id,
            fresh=fresh,
            use_ssd=use_ssd,
            scratch_dir=scratch_dir,
        )
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
    run_label = "Force Run Fresh" if fresh else "Force Run"
    return OrderSubmitResponse(
        status="accepted",
        order_id=order_id,
        service_code=job.service_code,
        message=f"Queued ({run_label} — {plugin.display_name if plugin else job.service_code})",
        queue_position=queue_position,
    )


_REPROCESS_SUPPORTED = frozenset({"carrier_screening", "whole_exome", "health_screening", "sgnipt"})


@app.post("/order/{order_id}/reprocess-results")
async def reprocess_order_results(order_id: str):
    """
    파이프라인 없이 process_results 만 다시 실행합니다 (Portal «Reprocess only»).

    - carrier_screening: 기존 VCF → annotation/QC/result.json 재생성
    - sgnipt: 기존 파이프라인 결과 JSON → result.json 재생성

    전체를 FASTQ부터 돌리려면 POST /order/{id}/start («Force Run»)를 사용하세요.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    if job.service_code not in _REPROCESS_SUPPORTED:
        raise HTTPException(
            status_code=400,
            detail=f"reprocess-results is not supported for service: {job.service_code}",
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
        raise HTTPException(status_code=400, detail=f"No plugin for {job.service_code}")

    # 서비스별 경로 정규화
    if job.service_code in _CARRIER_LIKE:
        from .services.carrier_screening.layout_norm import apply_carrier_layout_directories
        if apply_carrier_layout_directories(job):
            await queue_manager.persist_job(job)
    elif job.service_code == "sgnipt":
        from .services.sgnipt import apply_sgnipt_layout_directories
        if apply_sgnipt_layout_directories(job):
            await queue_manager.persist_job(job)

    completion_ok = await plugin.check_completion(job)
    if not completion_ok:
        detail_map = {
            "carrier_screening": "No VCF found for this order — cannot reprocess (check analysis/output paths).",
            "whole_exome": "No VCF found for this order — cannot reprocess (check analysis/output paths).",
            "health_screening": "No VCF found for this order — cannot reprocess (check analysis/output paths).",
            "sgnipt": "No pipeline result JSON found for this order — cannot reprocess (check output path).",
        }
        raise HTTPException(
            status_code=400,
            detail=detail_map.get(job.service_code, "Completion check failed — no output found."),
        )

    logger.info("[reprocess-results] Starting process_results for %s (%s)", order_id, job.service_code)
    process_ok = await plugin.process_results(job)
    if not process_ok:
        hint = (getattr(job, "error_log", None) or "").strip()
        detail = (
            f"process_results failed: {hint}"
            if hint
            else "process_results failed — see daemon logs (container: docker logs service-daemon)"
        )
        raise HTTPException(status_code=500, detail=detail)

    await queue_manager.finalize_reprocess_results(job)
    logger.info("[reprocess-results] Done for %s", order_id)

    return {
        "status": "ok",
        "order_id": order_id,
        "message": "Reprocessed (result.json updated)",
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


@app.get("/api/portal/bam-csv/browse")
async def browse_bam_csv_portal(
    path: str = Query(default="", description="BAM 데이터 루트 아래 상대 경로"),
    service_code: Optional[str] = Query(default=None, description="sgnipt | carrier_screening"),
    abs_path: Optional[str] = Query(default=None, description="임의 절대 경로 직접 탐색"),
):
    """BAM samplesheet CSV 파일 브라우져 (포털용 — 인증 면제 경로)."""
    return _browse_bam_csv_directory_payload(path, service_code, abs_path)


@app.get("/api/bam-csv/browse")
async def browse_bam_csv(
    path: str = Query(default="", description="BAM 데이터 루트 아래 상대 경로"),
    service_code: Optional[str] = Query(default=None, description="sgnipt | carrier_screening"),
    abs_path: Optional[str] = Query(default=None, description="임의 절대 경로 직접 탐색"),
):
    """BAM samplesheet CSV 파일 브라우져."""
    return _browse_bam_csv_directory_payload(path, service_code, abs_path)


@app.get("/api/portal/bam-csv/sample-ids")
async def read_bam_csv_sample_ids(
    abs_path: str = Query(description="BAM samplesheet CSV 절대 경로"),
):
    """CSV 첫 번째 컬럼(sample_id)을 읽어 반환합니다 (포털용 — 인증 면제)."""
    path_clean = abs_path.strip()
    if not os.path.isfile(path_clean):
        raise HTTPException(status_code=404, detail=f"File not found: {path_clean}")
    sample_ids: List[str] = []
    try:
        import csv as _csv
        with open(path_clean, newline="", encoding="utf-8-sig") as fh:
            reader = _csv.DictReader(fh)
            for row in reader:
                sid = (row.get("sample_id") or "").strip()
                if sid:
                    sample_ids.append(sid)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Cannot parse CSV: {e}") from e
    return {"abs_path": path_clean, "sample_ids": sample_ids}


@app.get("/api/portal/wes-panels")
async def portal_wes_panels():
    """
    Orderable WES/exome interpretation panels (disease BED + metadata).
    Bundled: ``data/wes_panels.json``. Custom (portal): ``data/wes_panels_custom.json`` — merged; custom overrides same id.
    """
    from app.services.wes_panels import load_wes_panels_raw, panels_for_api_response

    raw = load_wes_panels_raw()
    ver = raw.get("version", 1) if isinstance(raw, dict) else 1
    return {"version": ver, "panels": panels_for_api_response()}


@app.post("/api/portal/wes-panels/custom")
async def portal_wes_panels_custom_save(body: WesPanelCustomSave):
    """Create or update a user-defined panel package (saved under ``data/wes_panels_custom.json``)."""
    from app.services.wes_panels import save_custom_panel

    try:
        extra = save_custom_panel(
            body.id,
            body.label,
            body.category,
            body.description,
            body.backbone_bed,
            body.disease_bed,
            body.genes,
            body.genes_text,
            body.gene_source_bed,
            body.skip_generated_bed,
        )
        return {"status": "ok", **extra}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except OSError as e:
        raise HTTPException(
            status_code=503,
            detail=f"Could not write panel catalog file (check disk permissions and WES_PANELS_CUSTOM_JSON path): {e}",
        ) from e


@app.delete("/api/portal/wes-panels/custom/{panel_id:path}")
async def portal_wes_panels_custom_delete(panel_id: str):
    """Remove a custom panel (built-in ids cannot be deleted)."""
    from app.services.wes_panels import delete_custom_panel

    pid = (panel_id or "").strip().strip("/")
    if not pid:
        raise HTTPException(status_code=400, detail="panel_id required")
    try:
        ok = delete_custom_panel(pid)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    if not ok:
        raise HTTPException(status_code=404, detail="Custom panel not found")
    return {"status": "ok", "deleted": pid}


@app.get("/api/wes-panels")
async def api_wes_panels_alias():
    """Same payload as ``/api/portal/wes-panels`` (for clients that omit the portal prefix)."""
    return await portal_wes_panels()


@app.get("/api/portal/resources")
async def get_system_resources():
    """시스템 리소스 현황 (CPU/Memory/Disk) — 포털 대시보드용 (인증 면제)."""
    try:
        import psutil  # type: ignore
    except ImportError:
        raise HTTPException(status_code=503, detail="psutil not installed — run: pip install psutil")

    # CPU per-core (non-blocking: interval=None uses accumulated since last call)
    cpu_per_core = psutil.cpu_percent(interval=None, percpu=True)
    cpu_total = psutil.cpu_percent(interval=None)
    cpu_freq = psutil.cpu_freq(percpu=False)
    cpu_count_logical = psutil.cpu_count(logical=True)
    cpu_count_physical = psutil.cpu_count(logical=False)

    # Memory
    vm = psutil.virtual_memory()
    swap = psutil.swap_memory()

    # Disk — report each unique mount point visible inside the container
    _seen: set = set()
    disk_list = []
    priority_paths = ["/", "/data", "/home", "/tmp"]
    candidate_paths = priority_paths + [p.mountpoint for p in psutil.disk_partitions(all=False)]
    for mp in candidate_paths:
        if not os.path.isdir(mp):
            continue
        try:
            usage = psutil.disk_usage(mp)
        except (PermissionError, OSError):
            continue
        key = (usage.total, usage.used)
        if key in _seen:
            continue
        _seen.add(key)
        disk_list.append({
            "mount": mp,
            "total": usage.total,
            "used": usage.used,
            "free": usage.free,
            "percent": usage.percent,
        })

    return {
        "cpu": {
            "total_percent": cpu_total,
            "per_core": cpu_per_core,
            "logical_cores": cpu_count_logical,
            "physical_cores": cpu_count_physical,
            "freq_mhz": round(cpu_freq.current, 1) if cpu_freq else None,
        },
        "memory": {
            "total": vm.total,
            "available": vm.available,
            "used": vm.used,
            "percent": vm.percent,
        },
        "swap": {
            "total": swap.total,
            "used": swap.used,
            "free": swap.free,
            "percent": swap.percent,
        },
        "disk": disk_list,
    }


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

    if job.service_code in _CARRIER_LIKE:
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
                    "PDF report is not supported for this order configuration. "
                    "Set package / interpretation panel so the order resolves to carrier (standard/couples), "
                    "whole exome, proactive health, or PGx (e.g. WES panel category pgx / proactive_health, "
                    "or package_code PGx / WholeExome / HealthScreening)."
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
    try:
        success = await plugin.generate_report(
            job=job,
            confirmed_variants=request.confirmed_variants,
            reviewer_info=request.reviewer_info,
            patient_info=request.patient_info,
            partner_info=request.partner_info,
            languages=request.languages,
        )
    except Exception as e:
        logger.exception("Report generation raised for order %s", order_id)
        msg = str(e).strip() or repr(e)
        raise HTTPException(
            status_code=500,
            detail=msg[:8000],
        ) from e

    if not success:
        raise HTTPException(
            status_code=500,
            detail="Report generation failed (plugin returned false — check daemon logs).",
        )

    await queue_manager.mark_report_ready(
        order_id,
        message="Report ready for download",
    )

    store = queue_manager.store
    if store:
        review_langs = request.languages or []
        if job.service_code in _CARRIER_LIKE:
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
    if job.service_code in _CARRIER_LIKE:
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


# ─── Report Preview & From-HTML Endpoints ─────────────────

@app.post("/order/{order_id}/report/preview")
async def preview_report_html(order_id: str, request: ReportGenerateRequest):
    """
    Render the report Jinja template as HTML and return it (no PDF, no disk writes).
    Portal uses this for in-browser preview + manual editing before final generation.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    plugin = get_plugin(job.service_code)
    if not plugin or not hasattr(plugin, "generate_report"):
        raise HTTPException(status_code=400, detail="Service does not support report generation")

    if job.service_code not in _CARRIER_LIKE:
        raise HTTPException(status_code=400, detail="Preview is only supported for carrier-like services")

    from .services.carrier_screening.report import (
        carrier_report_template_kind,
        report_languages_from_order,
        _carrier_order_flat,
        generate_report_json,
        generate_report_pdf,
        _render_html_for_language,
        carrier_pdf_jinja_stem,
    )
    from .services.carrier_screening.plugin import (
        carrier_report_output_dir,
        resolve_carrier_pdf_template_dir,
        _extra_result_json_paths_for_carrier_report,
    )
    from .services.carrier_screening.report import CARRIER_PDF_SOLO_KINDS

    p_raw = job.params or {}
    kind = carrier_report_template_kind(p_raw)
    if kind is None:
        raise HTTPException(status_code=400, detail="PDF report is not supported for this order type.")
    langs_obj = report_languages_from_order(p_raw)
    if langs_obj is None:
        raise HTTPException(status_code=400, detail="Report language must be EN, CN, or KO.")
    languages = langs_obj if isinstance(langs_obj, list) else [langs_obj]

    output_dir = carrier_report_output_dir(job)
    template_dir = resolve_carrier_pdf_template_dir()
    params = _carrier_order_flat(p_raw)

    if kind in CARRIER_PDF_SOLO_KINDS:
        partner_info = None
    else:
        partner_info = request.partner_info

    pi = dict(request.patient_info) if request.patient_info else {}
    if not (pi.get("name") or "").strip():
        pi["name"] = (params.get("patient_name") or "").strip() or job.sample_name

    confirmed_variants = request.confirmed_variants or []
    report_json_path = os.path.join(output_dir, "report.json")

    def _render_preview():
        generate_report_json(
            order_id=job.order_id,
            sample_name=job.sample_name,
            confirmed_variants=confirmed_variants,
            reviewer_info=request.reviewer_info,
            qc_summary={},
            output_dir=output_dir,
            patient_info=pi,
            partner_info=partner_info,
            order_params=p_raw,
            report_language=(languages[0] if languages else "EN"),
            disease_gene_json=None,
            gene_knowledge_db=settings.gene_knowledge_db or None,
            gene_knowledge_enrich_on_report=False,
            gene_knowledge_gemini_on_report=False,
            gemini_api_key=None,
            gene_knowledge_gemini_model=None,
            extra_result_json_paths=_extra_result_json_paths_for_carrier_report(job),
            pdf_template_kind=kind,
        )

        with open(report_json_path, "r", encoding="utf-8") as f:
            report_data = json.load(f)

        from .services.carrier_screening.dark_genes import sanitize_dark_genes_payload_for_pdf_render
        from .services.carrier_screening.pgx_report import sanitize_pgx_payload_for_pdf_render
        from .services.carrier_screening.report import (
            _merge_dark_genes_from_result_json_for_pdf,
            _merge_pgx_from_result_json_for_pdf,
        )
        extras = _extra_result_json_paths_for_carrier_report(job)
        _merge_dark_genes_from_result_json_for_pdf(report_data, report_json_path, extra_result_json_paths=extras)
        sanitize_dark_genes_payload_for_pdf_render(report_data)
        _merge_pgx_from_result_json_for_pdf(report_data, report_json_path, extra_result_json_paths=extras)
        sanitize_pgx_payload_for_pdf_render(report_data)

        is_couple = report_data.get("report_metadata", {}).get("is_couple", False)
        result = {}
        for lang in languages:
            html_content = _render_html_for_language(report_data, lang, template_dir, is_couple)
            pdf_kind = report_data.get("report_metadata", {}).get("pdf_template_kind")
            stem = carrier_pdf_jinja_stem(pdf_kind, is_couple)
            result[lang] = {
                "html": html_content,
                "template": f"{stem}_{lang.upper()}.html",
            }
        return result

    try:
        rendered = await asyncio.to_thread(_render_preview)
    except Exception as e:
        logger.error(f"Report preview failed for {order_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Preview rendering failed: {e}")

    return JSONResponse({
        "status": "ok",
        "order_id": order_id,
        "languages": rendered,
    })


@app.post("/order/{order_id}/report/from-html")
async def generate_report_from_html(order_id: str, request: Request):
    """
    Accept edited HTML per language, generate PDFs, mark order REPORT_READY,
    and trigger platform upload — same finalization as POST /order/{id}/report.

    Body: { "languages": { "EN": "<html>...</html>", "CN": "<html>..." } }
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    if job.service_code not in _CARRIER_LIKE:
        raise HTTPException(status_code=400, detail="Only carrier-like services supported")

    body = await request.json()
    languages_html: dict = body.get("languages") or {}
    if not languages_html:
        raise HTTPException(status_code=400, detail="languages dict is required")

    from .services.carrier_screening.plugin import (
        carrier_report_output_dir,
        resolve_carrier_pdf_template_dir,
    )

    output_dir = carrier_report_output_dir(job)
    template_dir = resolve_carrier_pdf_template_dir()
    os.makedirs(output_dir, exist_ok=True)

    raw_name = (job.params or {}).get("patient_name") or job.sample_name or "Patient"
    patient_name = "".join([c if c.isalnum() else "_" for c in raw_name]).replace("__", "_")

    generated_pdfs = []

    def _render_all():
        from weasyprint import HTML as WeasyprintHTML
        base_url = template_dir if template_dir else output_dir
        for lang, html_content in languages_html.items():
            lang_u = (lang or "EN").strip().upper()
            html_content = (html_content or "").strip()
            if not html_content:
                continue
            html_filename = f"Report_{order_id}_{patient_name}_{lang_u}.html"
            pdf_filename = f"Report_{order_id}_{patient_name}_{lang_u}.pdf"
            html_path = os.path.join(output_dir, html_filename)
            pdf_path = os.path.join(output_dir, pdf_filename)
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_content)
            WeasyprintHTML(string=html_content, base_url=base_url).write_pdf(pdf_path)
            generated_pdfs.append(pdf_filename)
            logger.info(f"Generated PDF from edited HTML for {order_id}: {pdf_filename}")

    try:
        await asyncio.to_thread(_render_all)
    except ImportError:
        raise HTTPException(status_code=500, detail="weasyprint not installed")
    except Exception as e:
        logger.error(f"PDF from HTML failed for {order_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"PDF generation failed: {e}")

    if not generated_pdfs:
        raise HTTPException(status_code=400, detail="No PDFs were generated (empty HTML?)")

    await queue_manager.mark_report_ready(order_id, message="Report ready for download")

    store = queue_manager.store
    if store and output_dir:
        await asyncio.to_thread(ingest_report_json_from_disk, store, job, output_dir)

    report_files: list = []
    report_json_path = os.path.join(output_dir, "report.json")
    if os.path.exists(report_json_path):
        report_files.append(OutputFile(
            file_path=report_json_path, file_type="report_json",
            file_name="report.json", content_type="application/json",
        ))
    for pat, ftype, mime in [
        ("Report_*.pdf", "report_pdf", "application/pdf"),
        ("Report_*.html", "report_html", "text/html"),
    ]:
        for fp in glob.glob(os.path.join(output_dir, pat)):
            report_files.append(OutputFile(
                file_path=fp, file_type=ftype,
                file_name=os.path.basename(fp), content_type=mime,
            ))

    platform_client = get_platform_client()
    if report_files and settings.platform_api_enabled:
        async def _upload_bg():
            try:
                results = await platform_client.upload_all_outputs(
                    order_id, job.service_code, report_files
                )
                ok = sum(1 for r in results.values() if r.status.value == "SUCCESS")
                logger.info(f"Platform upload for {order_id}: {ok}/{len(report_files)}")
            except Exception as exc:
                logger.exception(f"Platform upload failed for {order_id}: {exc}")
        asyncio.create_task(_upload_bg())

    return JSONResponse({
        "status": "ok",
        "order_id": order_id,
        "report_files": generated_pdfs,
        "message": f"Report generated: {len(generated_pdfs)} PDF(s). Order marked REPORT_READY.",
    })


# ─── Order List & Result Endpoints ────────────────────────

def _job_to_order_summary_dict(j) -> Dict[str, Any]:
    """Same shape as each row in GET /orders (Portal list + follow-up copy)."""
    return {
        "order_id": j.order_id,
        "service_code": j.service_code,
        "sample_name": j.sample_name,
        "work_dir": j.work_dir,
        "status": j.status.value,
        "progress": j.progress,
        "message": j.message,
        "error_log": j.error_log,
        "created_at": j.created_at,
        "started_at": j.started_at,
        "completed_at": j.completed_at,
        "updated_at": getattr(j, "updated_at", None),
        "params": j.params or {},
        "fastq_r1_url": j.fastq_r1_url,
        "fastq_r2_url": j.fastq_r2_url,
        "fastq_r1_path": j.fastq_r1_path,
        "fastq_r2_path": j.fastq_r2_path,
    }


@app.get("/order/{order_id}")
async def get_order(order_id: str):
    """
    단일 주문 스냅샷 (목록 행과 동일 필드).

    Portal «New order from this» 등에서 전체 params 를 안전하게 복사할 때 사용합니다.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    return _job_to_order_summary_dict(job)


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
        "orders": [_job_to_order_summary_dict(j) for j in all_jobs],
        "total": len(all_jobs),
    }


@app.get("/order/{order_id}/result")
async def get_order_result(
    order_id: str,
    debug: bool = Query(False, description="Include _review_debug with disk paths (portal Review)."),
):
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

    # 경로가 DB에 저장되지 않은 경우 현재 설정으로 보정
    if job.service_code == "sgnipt" and not job.output_dir:
        from app.services.sgnipt import apply_sgnipt_layout_directories
        if apply_sgnipt_layout_directories(job):
            await queue_manager.persist_job(job)
    elif job.service_code in _CARRIER_LIKE and not job.output_dir:
        from app.services.carrier_screening.layout_norm import apply_carrier_layout_directories
        if apply_carrier_layout_directories(job):
            await queue_manager.persist_job(job)

    result_json_path = None
    if job.service_code in _CARRIER_LIKE:
        from app.services.carrier_screening.plugin import carrier_result_json_path

        result_json_path = carrier_result_json_path(job)
    elif job.output_dir:
        result_json_path = os.path.join(job.output_dir, "result.json")

    if result_json_path and os.path.isfile(result_json_path):
        with open(result_json_path, "r", encoding="utf-8") as f:
            result_data = json.load(f)
        # sgNIPT: process_results() 가 이미 variant_report 를 병합한 완성본을 저장합니다.
        # 구버전 result.json (clinical_findings 없는 요약 전용) 에 대한 즉석 fallback 병합.
        if isinstance(result_data, dict) and job.service_code == "sgnipt" \
                and not result_data.get("clinical_findings") \
                and result_data.get("samples"):
            logger.info("[sgnipt] result.json has no clinical_findings — running on-the-fly merge (legacy file)")
            from app.services.sgnipt import SgNIPTPlugin
            samples = result_data.get("samples") or []
            sample0 = samples[0] if samples else {}
            sample_id = sample0.get("sample_id") or (job.order_id or "").strip()
            vr_path = SgNIPTPlugin._find_variant_report_json(job.analysis_dir, job.output_dir, sample_id)
            vr: dict = {}
            if vr_path:
                try:
                    with open(vr_path, "r", encoding="utf-8") as _f:
                        vr = json.load(_f) or {}
                except Exception as _e:
                    logger.warning("[sgnipt] Could not load variant_report.json %s: %s", vr_path, _e)
            result_data = SgNIPTPlugin._build_merged_result(result_data, vr)
        if isinstance(result_data, dict):
            from app.services import get_plugin
            from app.services.carrier_screening.dark_genes import (
                ensure_dark_genes_detailed_sections,
                merge_visual_evidence_across_roots,
            )

            result_data = ensure_dark_genes_detailed_sections(dict(result_data))
            if job.service_code in _CARRIER_LIKE:
                pl = get_plugin(job.service_code)
                dg = result_data.get("dark_genes")
                if (
                    pl is not None
                    and hasattr(pl, "dark_genes_search_roots")
                    and isinstance(dg, dict)
                    and dg.get("status") not in ("not_found", "error")
                ):
                    roots = pl.dark_genes_search_roots(job)
                    dg2 = dict(dg)
                    dg2["visual_evidence"] = merge_visual_evidence_across_roots(roots)
                    result_data = {**result_data, "dark_genes": dg2}
            if job.service_code in _CARRIER_LIKE:
                from app.services.carrier_screening.review import (
                    enrich_review_build_with_smaca_vcf_depths,
                )

                meta = result_data.get("metadata")
                if isinstance(meta, dict):
                    rb = meta.get("review_build")
                    if isinstance(rb, dict) and not rb.get("smaca_snp_depths"):
                        enrich_review_build_with_smaca_vcf_depths(
                            rb, str(getattr(job, "sample_name", "") or "")
                        )
        if store and isinstance(result_data, dict):
            try:
                await asyncio.to_thread(store.set_result_json, order_id, result_data)
            except Exception as _store_err:
                logger.warning("Could not cache result_json for %s: %s", order_id, _store_err)
        if isinstance(result_data, dict):
            from app.services.carrier_screening.review import normalize_variants_for_portal

            payload = dict(result_data)
            payload["order_params"] = job.params or {}
            payload["service_code"] = job.service_code
            if isinstance(payload.get("variants"), list):
                payload["variants"] = normalize_variants_for_portal(payload["variants"])
            if debug:
                payload["_review_debug"] = _order_result_review_debug(
                    job, result_json_path, "disk", result_data
                )
            return JSONResponse(content=payload, headers=_RESULT_JSON_CACHE_HEADERS)
        extra = {"data": result_data, "order_params": job.params or {}}
        if debug:
            extra["_review_debug"] = _order_result_review_debug(job, result_json_path, "disk", None)
        return JSONResponse(
            content=extra,
            headers=_RESULT_JSON_CACHE_HEADERS,
        )

    if store:
        cached = await asyncio.to_thread(store.get_result_json, order_id)
        # Carrier: SQLite row can hold a snapshot copied from another order's result.json
        # (embedded order_id in the JSON body does not match this row's order_id).
        if (
            cached is not None
            and job.service_code in _CARRIER_LIKE
            and isinstance(cached, dict)
        ):
            body_oid = (cached.get("order_id") or "").strip()
            if body_oid and body_oid != (order_id or "").strip():
                logger.warning(
                    "Discarding stale result_json snapshot for order %s: "
                    "embedded order_id=%r does not match (likely cached from wrong file path).",
                    order_id,
                    body_oid,
                )
                try:
                    await asyncio.to_thread(store.clear_result_json, order_id)
                except Exception as _clr:
                    logger.warning("Could not clear stale result_json for %s: %s", order_id, _clr)
                cached = None
        if cached is not None:
            if isinstance(cached, dict):
                from app.services import get_plugin
                from app.services.carrier_screening.dark_genes import (
                    ensure_dark_genes_detailed_sections,
                    merge_visual_evidence_across_roots,
                )

                out = ensure_dark_genes_detailed_sections(dict(cached))
                if job.service_code in _CARRIER_LIKE:
                    pl = get_plugin(job.service_code)
                    dg = out.get("dark_genes")
                    if (
                        pl is not None
                        and hasattr(pl, "dark_genes_search_roots")
                        and isinstance(dg, dict)
                        and dg.get("status") not in ("not_found", "error")
                    ):
                        roots = pl.dark_genes_search_roots(job)
                        dg2 = dict(dg)
                        dg2["visual_evidence"] = merge_visual_evidence_across_roots(roots)
                        out = {**out, "dark_genes": dg2}
                if job.service_code in _CARRIER_LIKE:
                    from app.services.carrier_screening.review import (
                        enrich_review_build_with_smaca_vcf_depths,
                    )

                    meta = out.get("metadata")
                    if isinstance(meta, dict):
                        rb = meta.get("review_build")
                        if isinstance(rb, dict) and not rb.get("smaca_snp_depths"):
                            enrich_review_build_with_smaca_vcf_depths(
                                rb, str(getattr(job, "sample_name", "") or "")
                            )
                from app.services.carrier_screening.review import normalize_variants_for_portal

                payload = dict(out)
                payload["order_params"] = job.params or {}
                if isinstance(payload.get("variants"), list):
                    payload["variants"] = normalize_variants_for_portal(payload["variants"])
                if debug:
                    rp_path = None
                    if job.service_code in _CARRIER_LIKE:
                        from app.services.carrier_screening.plugin import carrier_result_json_path

                        rp_path = carrier_result_json_path(job)
                    elif job.output_dir:
                        rp_path = os.path.join(job.output_dir, "result.json")
                    payload["_review_debug"] = _order_result_review_debug(
                        job, rp_path, "db_snapshot", out
                    )
                return JSONResponse(content=payload, headers=_RESULT_JSON_CACHE_HEADERS)
            extra = {"data": cached, "order_params": job.params or {}}
            if debug:
                rp_path = None
                if job.service_code in _CARRIER_LIKE:
                    from app.services.carrier_screening.plugin import carrier_result_json_path

                    rp_path = carrier_result_json_path(job)
                elif job.output_dir:
                    rp_path = os.path.join(job.output_dir, "result.json")
                extra["_review_debug"] = _order_result_review_debug(
                    job, rp_path, "db_snapshot", None
                )
            return JSONResponse(
                content=extra,
                headers=_RESULT_JSON_CACHE_HEADERS,
            )

    detail = (
        f"result.json not found for order: {order_id} "
        f"(no file on disk for this order's review path and no valid DB snapshot)."
    )
    if job.service_code in _CARRIER_LIKE:
        try:
            from app.services.carrier_screening.plugin import carrier_report_output_dir

            detail += (
                f" Expected file: {os.path.join(carrier_report_output_dir(job), 'result.json')}"
                f" — run «Reprocess» or process_results after deploying per-order paths."
            )
        except Exception:
            pass
    raise HTTPException(status_code=404, detail=detail)


@app.get("/order/{order_id}/gene-coverage/{gene_symbol}")
async def get_order_gene_coverage(order_id: str, gene_symbol: str):
    """
    Twist exome target intervals for a gene (HGNC in BED column 4) plus optional per-gene depth
    sidecars — supports Variant Review (carrier / whole_exome / health_screening). Clinical panel
    is gene-list filtered; carrier disease/ACMG BEDs are not used for these intervals.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    if job.service_code not in _CARRIER_LIKE:
        raise HTTPException(
            status_code=400,
            detail=(
                "Gene panel coverage is only available for carrier_screening, whole_exome, "
                "and health_screening orders."
            ),
        )
    from app.services.gene_panel_coverage import build_gene_panel_coverage_report

    try:
        report = build_gene_panel_coverage_report(job, gene_symbol)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    return JSONResponse(content=report)


@app.get("/order/{order_id}/coverage-context")
async def get_order_coverage_context(order_id: str):
    """
    Portal **Coverage** tab: interpretation gene symbols (panel + extras) and BAM paths for IGV.js.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    _coverage_ok = _CARRIER_LIKE | {"sgnipt"}
    if job.service_code not in _coverage_ok:
        raise HTTPException(
            status_code=400,
            detail="Coverage context is only available for carrier_screening, whole_exome, health_screening, and sgnipt.",
        )
    if job.service_code == "sgnipt":
        genes: list = []
    else:
        from app.services.wes_panels import interpretation_gene_set_for_job
        genes = sorted(interpretation_gene_set_for_job(job))
    bams_all = _list_order_bam_tracks(job)
    # Never offer paraphase / SMN / FMR1 / … as auto-IGV tracks — only genome-wide candidates.
    bams_igv = [t for t in bams_all if not t.get("ancillary")]
    ctx: Dict[str, Any] = {
        "order_id": order_id,
        "interpretation_genes": genes,
        "bam_tracks": bams_igv,
        "genome_id": "hg38",
    }
    if bams_all and not bams_igv:
        ctx["igv_bam_message"] = (
            "Only ancillary BAMs (e.g. paraphase) were found under this order; they are not used for auto-IGV. "
            "Add the main exome BAM and .bai under analysis/output (or prior-reuse paths), or set "
            "params.igv_bam to that BAM’s absolute path."
        )
    if (job.params or {}).get("_prior_reuse"):
        pid = (job.params or {}).get("_prior_reuse_order_id")
        if isinstance(pid, str) and pid.strip():
            ctx["prior_reuse_order_id"] = pid.strip()
    return JSONResponse(
        content=ctx,
        headers=_RESULT_JSON_CACHE_HEADERS,
    )


async def _dark_genes_review_impl(order_id: str, body: DarkGenesReviewRequest) -> Dict[str, Any]:
    """
    Save per-section **Approve**, **Notes**, and **Risk** (PDF title accent) for the Dark genes detailed report.

    Writes ``dark_genes.section_reviews`` (index-aligned with ``detailed_sections``) into
    the same ``result.json`` that Generate Report reads (``carrier_report_output_dir``, mirroring
    to ``job.output_dir`` when those paths differ). Notes and workflow state are stored for the
    portal; the customer PDF lists all eligible supplementary sections (Overview / QC-only /
    duplicate row still omitted) and may show **Notes** when present.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    if job.service_code not in _CARRIER_LIKE:
        raise HTTPException(
            status_code=400,
            detail="dark-genes-review is only supported for carrier_screening, whole_exome, and health_screening",
        )
    from app.services.carrier_screening.plugin import (
        carrier_result_json_path,
        write_carrier_result_json_sync,
    )

    path = carrier_result_json_path(job)
    if not path:
        raise HTTPException(
            status_code=404,
            detail="result.json not found — run analysis / reprocess first",
        )

    def _apply() -> Dict[str, Any]:
        from app.services.carrier_screening.dark_genes import (
            align_section_reviews,
            ensure_dark_genes_detailed_sections,
            _coerce_risk_level,
        )

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        data = ensure_dark_genes_detailed_sections(data)
        dg = data.get("dark_genes")
        if not isinstance(dg, dict):
            raise ValueError("result.json has no dark_genes object")
        sections = dg.get("detailed_sections")
        if not isinstance(sections, list) or len(sections) == 0:
            raise ValueError(
                "No detailed_sections (and no detailed_text to parse) — reprocess results "
                "so dark_genes has supplementary detailed content."
            )
        n = len(sections)
        incoming = []
        for x in body.section_reviews:
            item = {
                "approved": bool(x.approved),
                "notes": (x.notes or "")[:8000],
            }
            if x.risk is not None:
                item["risk"] = _coerce_risk_level(x.risk)
            incoming.append(item)
        dg2 = dict(dg)
        dg2["section_reviews"] = align_section_reviews(incoming, n, sections)
        data["dark_genes"] = dg2
        write_carrier_result_json_sync(job, data)
        return data

    try:
        data = await asyncio.to_thread(_apply)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e

    store = queue_manager.store
    if store and isinstance(data, dict):
        await asyncio.to_thread(store.set_result_json, order_id, data)

    dg_out = data.get("dark_genes") if isinstance(data, dict) else {}
    return {
        "status": "ok",
        "order_id": order_id,
        "dark_genes": dg_out if isinstance(dg_out, dict) else {},
    }


@app.patch("/order/{order_id}/dark-genes-review")
async def patch_dark_genes_review(order_id: str, body: DarkGenesReviewRequest):
    return await _dark_genes_review_impl(order_id, body)


@app.post("/order/{order_id}/dark-genes-review")
async def post_dark_genes_review(order_id: str, body: DarkGenesReviewRequest):
    """Same as PATCH — POST is supported for proxies that block PATCH."""
    return await _dark_genes_review_impl(order_id, body)


async def _pgx_review_impl(order_id: str, body: PgxReviewRequest) -> Dict[str, Any]:
    """
    Save ``pgx.portal_review`` (reviewer notes, reviewed flag) into the same ``result.json``
    as Generate Report (``carrier_report_output_dir``, mirrored to ``job.output_dir``).
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    if job.service_code not in _CARRIER_LIKE:
        raise HTTPException(
            status_code=400,
            detail="pgx-review is only supported for carrier_screening, whole_exome, and health_screening",
        )
    from app.services.carrier_screening.plugin import (
        carrier_result_json_path,
        write_carrier_result_json_sync,
    )

    path = carrier_result_json_path(job)
    if not path:
        raise HTTPException(
            status_code=404,
            detail="result.json not found — run analysis / reprocess first",
        )

    def _apply() -> Dict[str, Any]:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            data = {}
        pgx = data.get("pgx")
        if not isinstance(pgx, dict):
            pgx = {"status": "not_found", "message": "PGx block was missing; created by portal save"}
        pgx2 = dict(pgx)
        prev_pr = pgx2.get("portal_review") if isinstance(pgx2.get("portal_review"), dict) else {}
        pgx2["portal_review"] = {
            **prev_pr,
            "reviewer_notes": (body.reviewer_notes or "")[:16000],
            "reviewed": bool(body.reviewed),
            "include_apoe_proactive_pdf": bool(body.include_apoe_proactive_pdf),
        }
        grs = pgx2.get("gene_results")
        if isinstance(grs, list) and body.gene_reviews is not None:
            by_gene = {x.gene: x for x in body.gene_reviews}
            new_grs = []
            for row in grs:
                if not isinstance(row, dict):
                    continue
                r = dict(row)
                g = r.get("gene")
                if g in by_gene:
                    u = by_gene[g]
                    r["reviewer_confirmed"] = bool(u.reviewer_confirmed)
                    r["reviewer_comment"] = (u.reviewer_comment or "")[:4000]
                new_grs.append(r)
            pgx2["gene_results"] = new_grs
        cgrs = pgx2.get("custom_gene_results")
        if isinstance(cgrs, list) and body.custom_gene_reviews is not None:
            by_key = {(x.gene, x.rsid): x for x in body.custom_gene_reviews}
            new_cgrs = []
            for row in cgrs:
                if not isinstance(row, dict):
                    continue
                r = dict(row)
                key = (r.get("gene", ""), r.get("rsid", ""))
                if key in by_key:
                    u = by_key[key]
                    r["reviewer_confirmed"] = bool(u.reviewer_confirmed)
                    r["reviewer_comment"] = (u.reviewer_comment or "")[:4000]
                new_cgrs.append(r)
            pgx2["custom_gene_results"] = new_cgrs
        data["pgx"] = pgx2
        write_carrier_result_json_sync(job, data)
        return data

    data = await asyncio.to_thread(_apply)

    store = queue_manager.store
    if store and isinstance(data, dict):
        await asyncio.to_thread(store.set_result_json, order_id, data)

    pgx_out = data.get("pgx") if isinstance(data, dict) else {}
    return {
        "status": "ok",
        "order_id": order_id,
        "pgx": pgx_out if isinstance(pgx_out, dict) else {},
    }


@app.patch("/order/{order_id}/pgx-review")
async def patch_pgx_review(order_id: str, body: PgxReviewRequest):
    return await _pgx_review_impl(order_id, body)


@app.post("/order/{order_id}/pgx-review")
async def post_pgx_review(order_id: str, body: PgxReviewRequest):
    """Same as PATCH — POST is supported for proxies that block PATCH."""
    return await _pgx_review_impl(order_id, body)


def _load_result_dict_for_gene_knowledge(order_id: str, job) -> Optional[Dict[str, Any]]:
    """Read result.json-like dict (variants list) without order_params merge."""
    result_json_path = None
    if getattr(job, "service_code", None) in _CARRIER_LIKE:
        from app.services.carrier_screening.plugin import carrier_result_json_path

        result_json_path = carrier_result_json_path(job)
    elif job.output_dir:
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


def _order_artifact_roots(job: Job) -> List[str]:
    """
    Directories to search for order output files (same order as download resolution).

    For carrier-like services, ``carrier_report_output_dir`` is listed first when set so
    ``Report_*.pdf`` / ``report.json`` match Generate Report even when
    ``CARRIER_SCREENING_REPORT_OUTPUT_ROOT`` differs from ``job.output_dir``.

    When ``_prior_reuse`` is set (reflex / new order from prior), the prior run's
    analysis/output trees are prepended so BAMs, mosdepth, and ``/order/.../file`` match
    the source sequencing run.
    """
    roots: List[str] = []
    roots.extend(prior_reuse_artifact_roots(job))
    if (job.service_code or "").strip() in _CARRIER_LIKE:
        from .services.carrier_screening.plugin import carrier_report_output_dir

        try:
            cro = carrier_report_output_dir(job)
            if cro:
                roots.append(cro)
        except Exception:
            pass
        # Nextflow outdir (IGV html, repeat SVGs under snapshots/ / repeat/)
        ad = (job.analysis_dir or "").strip()
        if ad:
            try:
                if os.path.isdir(ad):
                    roots.append(ad)
            except OSError:
                pass
    if job.output_dir:
        roots.append(job.output_dir)
    if (job.service_code or "").strip() in _CARRIER_LIKE:
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
    out: List[str] = []
    for root in roots:
        if not root:
            continue
        try:
            key = os.path.realpath(root)
        except OSError:
            continue
        if not os.path.isdir(root):
            continue
        if key in seen:
            continue
        seen.add(key)
        out.append(root)
    # Dirs used for mosdepth / gene coverage (plugin dark_genes roots, etc.) — needed so
    # BAM next to *.per-base.bed.gz resolves for /order/.../file and IGV.
    try:
        from .services.gene_panel_coverage import _coverage_search_roots

        for root in _coverage_search_roots(job):
            if not root:
                continue
            try:
                key = os.path.realpath(root)
            except OSError:
                continue
            if not os.path.isdir(root):
                continue
            if key in seen:
                continue
            seen.add(key)
            out.append(root)
    except Exception:
        pass
    return out


def _order_file_basename_from_path_or_url(path_or_url: Optional[str]) -> str:
    """Bare filename from a local path or http(s) URL path component."""
    if not path_or_url or not isinstance(path_or_url, str):
        return ""
    s = path_or_url.strip()
    if s.lower().startswith(("http://", "https://")):
        u = urlparse(s)
        return unquote(os.path.basename((u.path or "").rstrip("/")))
    return os.path.basename(s.replace("\\", "/"))


def _fastq_basename_alignment_match_stems(name: str) -> List[str]:
    """
    Substrings likely shared with the primary BAM: strip FASTQ extensions, drop Illumina
    read tokens, and add the prefix before _R1_/_R2_ (e.g. ..._S106 vs ..._R1_001_...).
    """
    out: List[str] = []
    base = (name or "").strip()
    if not base:
        return out
    low = base.lower()
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if low.endswith(suf):
            base = base[: -len(suf)].strip()
            low = base.lower()
            break

    def add(stem: str) -> None:
        stem = stem.strip().strip("._")
        if len(stem) >= 12:
            out.append(stem)

    ubase = base.upper()
    for marker in ("_R1_", "_R2_"):
        ix = ubase.find(marker.upper())
        if ix >= 12:
            add(base[:ix])

    stripped = re.sub(r"_R[12]_\d{3}(?=_|$)", "", base, flags=re.IGNORECASE).strip()
    add(stripped)
    slo = stripped.lower()
    j = slo.find("_downsampled")
    if j >= 12:
        add(stripped[:j])

    seen: set = set()
    uniq: List[str] = []
    for s in sorted(out, key=len, reverse=True):
        k = s.lower()
        if k not in seen:
            seen.add(k)
            uniq.append(s)
    return uniq


def _fastq_derived_bam_match_stems(job: Job) -> List[str]:
    stems: List[str] = []
    for raw in (job.fastq_r1_path, job.fastq_r2_path, job.fastq_r1_url, job.fastq_r2_url):
        bn = _order_file_basename_from_path_or_url(raw)
        stems.extend(_fastq_basename_alignment_match_stems(bn))
    seen: set = set()
    out: List[str] = []
    for s in stems:
        k = s.lower()
        if k not in seen:
            seen.add(k)
            out.append(s)
    out.sort(key=len, reverse=True)
    return out


def _order_bam_scan_path_ok(path: str) -> bool:
    norm = path.replace("\\", "/").lower()
    return "/env/" not in norm and "/viz_env/" not in norm


def _order_bam_search_roots(job: Job) -> List[str]:
    """
    Roots to scan for BAMs. Prefer pipeline ``analysis_dir`` / ``output_dir`` before
    ``carrier_report_output_dir`` so the primary exome/carrier alignments win over
    incidental copies (e.g. ``no_dup.bam``) under the report mirror tree.
    """
    seen: set = set()
    out: List[str] = []

    def add(path: Optional[str]) -> None:
        p = (path or "").strip()
        if not p or not os.path.isdir(p):
            return
        try:
            key = os.path.realpath(p)
        except OSError:
            return
        if key in seen:
            return
        seen.add(key)
        out.append(p)

    for p in prior_reuse_artifact_roots(job):
        add(p)
    add(job.analysis_dir)
    add(job.output_dir)
    for r in _order_artifact_roots(job):
        add(r)
    return out


_ANCILLARY_BAM_TOKENS = (
    "paraphase", "eh_realigned", "fmr1", "fragile", "frax",
    "smn_combined", "smn_merged", "smn_unified",
    "smn1_realigned", "smn2_realigned", "smn1_aln", "smn2_aln",
    "smn1_tagged", "smn2_tagged",
    "smn_raw", "strc_aln", "strc_realigned", "strc_gene2_realigned",
    "pms2_aln", "pms2_realigned",
    "gba_aln", "gba_realigned", "gba_tagged",
    "hba_aln", "hba_realigned",
    "no_dup", "nodup", "no-dup", "pre_dedup", "before_dedup",
)

_ANCILLARY_BAM_SUFFIXES = (
    "_realigned_old.bam", "_realigned_tagged.bam",
)


def _is_ancillary_bam(basename: str) -> bool:
    b = basename.lower()
    for tok in _ANCILLARY_BAM_TOKENS:
        if tok in b:
            return True
    for suf in _ANCILLARY_BAM_SUFFIXES:
        if b.endswith(suf):
            return True
    return False


def _bam_basename_preference_rank(job: Job, basename: str) -> int:
    """
    Lower = better candidate for IGV primary track.

    Ancillary / targeted BAMs (paraphase, SMN, FMR1, STRC, GBA, HBA, PMS2, no_dup, …)
    are checked **first** so a shared FASTQ stem cannot override the blacklist.
    """
    b = basename.lower()
    if _is_ancillary_bam(basename):
        return 8
    # Markdup / recal are the canonical pipeline output.
    if any(tok in b for tok in ("markdup", "mkdup", "recalibrated", "bqsr", ".md.")):
        return 0
    # Match main alignments named after the sample or order id.
    for token in (job.sample_name, job.order_id):
        t = (token or "").strip().lower().replace("-", "_")
        if t and len(t) >= 2 and t in b.replace("-", "_"):
            return 1
    # FASTQ library name on the order often matches the dedup/recal BAM stem.
    for stem in _fastq_derived_bam_match_stems(job):
        if stem.lower() in b:
            return 1
    return 3


def _bam_path_sort_key(job: Job, fp: str) -> Tuple[int, float]:
    rank = _bam_basename_preference_rank(job, os.path.basename(fp))
    try:
        mtime = os.path.getmtime(fp)
    except OSError:
        mtime = 0.0
    return (rank, -mtime)


def _safe_mtime(path: str) -> float:
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0.0


def _mosdepth_per_base_stem(per_base_basename: str) -> str:
    for suf in (
        ".mosdepth.global.per-base.bed.gz",
        ".mosdepth.per-base.bed.gz",
        ".global.per-base.bed.gz",
        ".per-base.bed.gz",
    ):
        if per_base_basename.endswith(suf):
            return per_base_basename[: -len(suf)]
    return ""


def _bam_paths_in_dir_sharing_mosdepth_stem(dir_path: str, stem: str) -> List[str]:
    if not stem or not os.path.isdir(dir_path):
        return []
    out_paths: List[str] = []
    sl = stem.lower()
    try:
        for fn in os.listdir(dir_path):
            low = fn.lower()
            if not low.endswith(".bam"):
                continue
            fp = os.path.join(dir_path, fn)
            if not os.path.isfile(fp):
                continue
            base = fn[: -4]
            bl = base.lower()
            if bl == sl or base.startswith(stem + ".") or base.startswith(stem + "_"):
                out_paths.append(fp)
    except OSError:
        return []
    seen_rp: set = set()
    uniq_paths: List[str] = []
    for p in out_paths:
        try:
            rp = os.path.realpath(p)
        except OSError:
            rp = p
        if rp in seen_rp:
            continue
        seen_rp.add(rp)
        uniq_paths.append(p)
    return uniq_paths


def _order_file_rel_from_abs_bam(job: Job, bam_abs: str) -> Optional[Dict[str, Any]]:
    """Build track dict for ``/order/{{id}}/file/{{rel}}`` if the BAM sits under a known root."""
    try:
        bam_abs = os.path.realpath(bam_abs)
    except OSError:
        return None
    if not os.path.isfile(bam_abs) or not bam_abs.lower().endswith(".bam"):
        return None
    for root in _order_artifact_roots(job):
        try:
            root_r = os.path.realpath(root)
        except OSError:
            continue
        if not (bam_abs == root_r or bam_abs.startswith(root_r + os.sep)):
            continue
        try:
            rel = os.path.relpath(bam_abs, root).replace("\\", "/")
        except ValueError:
            continue
        if ".." in rel:
            continue
        bai = bam_abs + ".bai"
        index_rel: Optional[str] = None
        if os.path.isfile(bai):
            try:
                index_rel = os.path.relpath(bai, root).replace("\\", "/")
            except ValueError:
                index_rel = None
        return {
            "rel_path": rel,
            "label": os.path.basename(bam_abs),
            "has_index": bool(index_rel),
            "index_rel_path": index_rel,
        }
    return None


def _pick_best_bam_path_for_igv(
    job: Job, paths: List[str], *, allow_ancillary_fallback: bool = False
) -> Optional[str]:
    """
    Pick one BAM path. By default, returns None if only ``rank >= 8`` (paraphase / SMN / …)
    candidates exist — caller should widen search instead of auto-loading an ancillary BAM.
    """
    if not paths:
        return None
    with_bai = [p for p in paths if os.path.isfile(p + ".bai")]
    pool = with_bai if with_bai else paths
    if not pool:
        return None
    non_anc = [p for p in pool if _bam_basename_preference_rank(job, os.path.basename(p)) < 8]
    use = non_anc if non_anc else (pool if allow_ancillary_fallback else [])
    if not use:
        return None
    use.sort(
        key=lambda p: (
            _bam_basename_preference_rank(job, os.path.basename(p)),
            -_safe_mtime(p),
        )
    )
    return use[0]


def _igv_bam_track_from_explicit_param(job: Job) -> Optional[Dict[str, Any]]:
    raw = (job.params or {}).get("igv_bam")
    if raw is None:
        car = (job.params or {}).get("carrier")
        if isinstance(car, dict):
            raw = car.get("igv_bam")
    if not isinstance(raw, str) or not raw.strip():
        return None
    p = os.path.abspath(os.path.normpath(raw.strip()))
    if not os.path.isfile(p):
        return None
    return _order_file_rel_from_abs_bam(job, p)


def _igv_bam_track_from_mosdepth_sibling(job: Job) -> Optional[Dict[str, Any]]:
    try:
        from .services.gene_panel_coverage import _coverage_search_roots, _find_mosdepth_per_base_bed
    except Exception:
        return None
    roots = _coverage_search_roots(job)
    pb = _find_mosdepth_per_base_bed(roots)
    if not pb or not os.path.isfile(pb):
        return None
    stem = _mosdepth_per_base_stem(os.path.basename(pb))
    if len(stem) < 2:
        return None
    d = os.path.dirname(pb)
    for _ in range(5):
        cands = _bam_paths_in_dir_sharing_mosdepth_stem(d, stem)
        best = _pick_best_bam_path_for_igv(job, cands, allow_ancillary_fallback=False)
        if best:
            return _order_file_rel_from_abs_bam(job, best)
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return None


def _list_order_bam_tracks_glob(job: Job, cap: int) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    seen_real: set = set()
    for root in _order_bam_search_roots(job):
        if len(out) >= cap:
            break
        if not root or not os.path.isdir(root):
            continue
        try:
            matches = glob.glob(os.path.join(root, "**", "*.bam"), recursive=True)
        except Exception:
            continue
        matches = [p for p in matches if os.path.isfile(p) and _order_bam_scan_path_ok(p)]
        try:
            matches.sort(key=lambda p: _bam_path_sort_key(job, p))
        except Exception:
            pass
        for fp in matches:
            if len(out) >= cap:
                break
            try:
                rk = os.path.realpath(fp)
            except OSError:
                continue
            if rk in seen_real:
                continue
            seen_real.add(rk)
            try:
                rel = os.path.relpath(fp, root).replace("\\", "/")
            except ValueError:
                continue
            if ".." in rel:
                continue
            bai = fp + ".bai"
            index_rel: Optional[str] = None
            if os.path.isfile(bai):
                try:
                    index_rel = os.path.relpath(bai, root).replace("\\", "/")
                except ValueError:
                    index_rel = None
            lab = os.path.basename(fp)
            rk = _bam_basename_preference_rank(job, lab)
            out.append(
                {
                    "rel_path": rel,
                    "label": lab,
                    "has_index": bool(index_rel),
                    "index_rel_path": index_rel,
                    "ancillary": rk >= 8,
                }
            )
    try:
        out.sort(
            key=lambda t: (
                1 if t.get("ancillary") else 0,
                _bam_basename_preference_rank(job, (t.get("label") or "")),
            )
        )
    except Exception:
        pass
    return out


def _annotate_bam_track_ancillary(job: Job, track: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not track:
        return None
    t = dict(track)
    lab = t.get("label") or ""
    t["ancillary"] = _bam_basename_preference_rank(job, lab) >= 8
    return t


def _list_order_bam_tracks(job: Job, cap: int = 32) -> List[Dict[str, Any]]:
    """
    BAM list for IGV: (1) optional ``params.igv_bam`` or ``params.carrier.igv_bam`` absolute path,
    (2) BAM in the same folder as indexed mosdepth ``*.per-base.bed.gz`` (matching file stem),
    (3) recursive scan with filename heuristics (FASTQ stem, markdup, deprioritize SMN/FMR1/paraphase).

    Each entry may include ``ancillary: true`` for paraphase / SMN / FMR1 / … — the portal skips these
    for auto-IGV when a non-ancillary indexed BAM exists.
    """
    primary: Optional[Dict[str, Any]] = None
    ex = _igv_bam_track_from_explicit_param(job)
    if ex and ex.get("has_index"):
        primary = ex
    else:
        md = _igv_bam_track_from_mosdepth_sibling(job)
        if md and md.get("has_index"):
            primary = md
    rest = _list_order_bam_tracks_glob(job, cap)
    if not primary:
        return rest
    primary = _annotate_bam_track_ancillary(job, primary)
    p_full = _resolve_order_artifact_path(job, str(primary.get("rel_path") or ""))
    try:
        p_key = os.path.realpath(p_full) if p_full else None
    except OSError:
        p_key = None
    if not p_key:
        return [primary] + rest[: cap - 1]
    filtered: List[Dict[str, Any]] = []
    for t in rest:
        rel = (t.get("rel_path") or "").strip()
        if not rel:
            continue
        other = _resolve_order_artifact_path(job, rel)
        try:
            ok = os.path.realpath(other) if other else None
        except OSError:
            ok = None
        if ok and ok == p_key:
            continue
        filtered.append(t)
    return [primary] + filtered[: cap - 1]


@app.get("/order/{order_id}/files")
async def list_order_files(order_id: str):
    """
    주문의 출력 파일 목록을 조회합니다.
    """
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")

    roots = _order_artifact_roots(job)
    if not roots:
        return JSONResponse(
            content={"files": [], "total": 0},
            headers={
                "Cache-Control": "no-store, no-cache, must-revalidate",
                "Pragma": "no-cache",
            },
        )

    # Same basename may exist under multiple roots (report output root vs job.output_dir).
    # Prefer the newest file by mtime so the PDF button matches a freshly generated report.
    best: Dict[str, Tuple[str, int, int]] = {}
    for output_dir in roots:
        try:
            names = os.listdir(output_dir)
        except OSError:
            continue
        for fname in names:
            fpath = os.path.join(output_dir, fname)
            if not os.path.isfile(fpath):
                continue
            try:
                mtime_ms = int(round(os.path.getmtime(fpath) * 1000))
            except OSError:
                mtime_ms = 0
            try:
                sz = os.path.getsize(fpath)
            except OSError:
                sz = 0
            prev = best.get(fname)
            if prev is None or mtime_ms > prev[1]:
                best[fname] = (fpath, mtime_ms, sz)

    files: List[Dict[str, Any]] = []
    for fname in sorted(best.keys()):
        fpath, mtime_ms, sz = best[fname]
        files.append({
            "name": fname,
            "size": sz,
            "mtime_ms": mtime_ms,
            "type": _guess_file_type(fname),
        })

    return JSONResponse(
        content={"files": files, "total": len(files)},
        headers={
            "Cache-Control": "no-store, no-cache, must-revalidate",
            "Pragma": "no-cache",
        },
    )


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
    Resolve a file under the same roots as ``list_order_files`` (see ``_order_artifact_roots``).

    If the same relative path exists in more than one root, return the path with the
    **newest** modification time (matches ``list_order_files`` mtime-based merge).
    """
    rel = (filename or "").strip().lstrip("/")
    if not rel or ".." in rel.replace("\\", "/"):
        return None
    hits: List[str] = []
    for root in _order_artifact_roots(job):
        hit = _safe_order_file_path(root, rel)
        if hit:
            hits.append(hit)
    if not hits:
        return None
    if len(hits) == 1:
        return hits[0]
    try:
        return max(hits, key=lambda p: os.path.getmtime(p))
    except OSError:
        return hits[0]


def _resolve_order_file_or_404(order_id: str, filename: str) -> str:
    queue_manager = get_queue_manager()
    job = queue_manager.get_job(order_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Order not found: {order_id}")
    file_path = _resolve_order_artifact_path(job, filename)
    if not file_path:
        raise HTTPException(status_code=404, detail=f"File not found: {filename}")
    return file_path


@app.head("/order/{order_id}/file/{filename:path}")
async def head_order_file(order_id: str, filename: str):
    """HEAD for IGV.js / byte-range clients that probe Content-Length before Range GETs."""
    file_path = _resolve_order_file_or_404(order_id, filename)
    try:
        stat = os.stat(file_path)
    except OSError:
        raise HTTPException(status_code=404, detail="File not found")
    return PlainTextResponse(
        content="",
        headers={
            "Content-Length": str(stat.st_size),
            "Content-Type": _guess_content_type(filename),
            "Accept-Ranges": "bytes",
            "Cache-Control": "no-store, no-cache, must-revalidate",
            "Pragma": "no-cache",
        },
    )


@app.get("/order/{order_id}/file/{filename:path}")
async def download_order_file(order_id: str, filename: str):
    """
    주문의 특정 출력 파일을 다운로드합니다.
    filename 에 슬래시 포함 가능 (예 qc/plot.png).
    """
    file_path = _resolve_order_file_or_404(order_id, filename)
    base_name = os.path.basename(filename) or "download"
    return FileResponse(
        path=file_path,
        filename=base_name,
        media_type=_guess_content_type(filename),
        headers={
            "Cache-Control": "no-store, no-cache, must-revalidate",
            "Pragma": "no-cache",
        },
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
        "svg": "image/svg+xml",
        "webp": "image/webp",
    }
    return ct_map.get(ext, "application/octet-stream")


# ─── Queue Endpoints ──────────────────────────────────────

@app.get("/queue/summary", response_model=QueueSummary)
async def get_queue_summary():
    """큐 상태 요약"""
    queue_manager = get_queue_manager()
    return queue_manager.get_summary()


def _dashboard_order_updated_iso(j: Job) -> str:
    """상태 갱신 시각(포털 Order updated): 우선 ``updated_at``."""
    return j.updated_at or j.completed_at or j.started_at or j.created_at or ""


@app.get("/queue/dashboard-bucket")
async def get_dashboard_bucket(
    bucket: str = Query(..., description="queued | running | completed | failed"),
    service_code: Optional[str] = Query(default=None, description="서비스 코드; 생략 시 전체"),
    sort: str = Query(default="order_updated", description="order_id | status | order_updated | message"),
    order: str = Query(default="desc", description="asc | desc"),
):
    """
    대시보드 통계 카드(버킷)에 대응하는 주문 목록.
    ``order_updated`` 는 Job 의 ``updated_at``(없으면 completed/started/created 로 대체)입니다.
    """
    allowed_bucket = {"queued", "running", "completed", "failed"}
    b = (bucket or "").strip().lower()
    if b not in allowed_bucket:
        raise HTTPException(status_code=400, detail=f"Invalid bucket: {bucket}")

    order_l = (order or "desc").strip().lower()
    if order_l not in ("asc", "desc"):
        order_l = "desc"
    sk = (sort or "order_updated").strip().lower()
    allowed_sort = {"order_id", "status", "order_updated", "message"}
    if sk not in allowed_sort:
        sk = "order_updated"

    queue_manager = get_queue_manager()
    jobs = queue_manager.get_dashboard_bucket_jobs(b, service_code)
    rev = order_l == "desc"

    if sk == "order_id":
        jobs = sorted(jobs, key=lambda j: (j.order_id or "").lower(), reverse=rev)
    elif sk == "status":
        jobs = sorted(
            jobs,
            key=lambda j: j.status.value if j.status else "",
            reverse=rev,
        )
    elif sk == "message":
        jobs = sorted(jobs, key=lambda j: (j.message or "").lower(), reverse=rev)
    else:
        jobs = sorted(jobs, key=_dashboard_order_updated_iso, reverse=rev)

    return {
        "bucket": b,
        "service_code": service_code,
        "sort": sk,
        "order": order_l,
        "total": len(jobs),
        "orders": [
            {
                "order_id": j.order_id,
                "status": j.status.value,
                "order_updated": _dashboard_order_updated_iso(j),
                "message": j.message or "",
            }
            for j in jobs
        ],
    }


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

def _accept_prefers_json(accept: str) -> bool:
    if not accept or not accept.strip():
        return False
    first = accept.split(",")[0].strip().split(";")[0].strip().lower()
    return first == "application/json"


@app.get("/")
async def dashboard(request: Request):
    """
    JSON dashboard for API clients and curl. Browsers (Accept: text/html first or present)
    are redirected to the portal unless JSON is explicitly preferred.
    """
    accept = request.headers.get("accept") or ""
    if not _accept_prefers_json(accept) and "text/html" in accept.lower():
        return RedirectResponse("/portal/", status_code=302)

    queue_manager = get_queue_manager()
    summary = queue_manager.get_summary()
    plugins = get_all_plugins()

    return {
        "service": "service-daemon",
        "version": "2.0.0",
        "environment": settings.app_env,
        "api_docs": {
            "swagger_ui": "/docs",
            "redoc": "/redoc",
            "openapi_json": "/openapi.json",
            "note": "Interactive HTML API reference (generated from OpenAPI). Prefer these over static markdown for exact paths and schemas.",
        },
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


# ─── Portal under /api/portal (proxies that only forward /api/* to this app) ───
_portal_dir_api = os.path.join(os.path.dirname(__file__), "..", "portal")
if os.path.isdir(_portal_dir_api):

    @app.get("/api/portal", include_in_schema=False)
    def portal_under_api_redirect():
        return RedirectResponse("/api/portal/")

    app.mount(
        "/api/portal",
        StaticFiles(directory=_portal_dir_api, html=True),
        name="portal_under_api",
    )
