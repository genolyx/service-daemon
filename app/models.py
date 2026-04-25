"""
Service Daemon Data Models

모든 서비스에 공통으로 사용되는 데이터 모델을 정의합니다.
기존 nipt-daemon의 models.py를 확장하여 service_code를 추가합니다.
"""

from enum import Enum
from typing import Optional, Dict, Any, List, Literal
from pydantic import BaseModel, Field, model_validator

from .datetime_kst import now_kst_iso


# ─── Enums ─────────────────────────────────────────────────

class OrderStatus(str, Enum):
    """주문 상태"""
    RECEIVED = "RECEIVED"           # 접수됨
    SAVED = "SAVED"                 # Portal에서만 저장됨, 파이프라인 미시작
    QUEUED = "QUEUED"               # 큐 대기 중
    DOWNLOADING = "DOWNLOADING"     # 입력 파일 다운로드 중
    RUNNING = "RUNNING"             # 파이프라인 실행 중
    PROCESSING = "PROCESSING"       # 결과 후처리 중 (annotation 등)
    UPLOADING = "UPLOADING"         # 결과 업로드 중
    COMPLETED = "COMPLETED"         # 분석 완료 (리뷰/리포트 전)
    REPORT_READY = "REPORT_READY"   # 최종 리포트 생성 완료 (Portal 다운로드 가능)
    FAILED = "FAILED"               # 실패
    CANCELLED = "CANCELLED"         # 취소됨


class ServiceCode(str, Enum):
    """지원 서비스 코드"""
    NIPT = "nipt"
    CARRIER_SCREENING = "carrier_screening"
    SG_NIPT = "sgnipt"
    WHOLE_EXOME = "whole_exome"
    HEALTH_SCREENING = "health_screening"


class NotificationStatus(str, Enum):
    """알림 상태"""
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    NOT_FOUND = "NOT_FOUND"
    SKIPPED = "SKIPPED"


# ─── Request / Response Models ─────────────────────────────

class OrderSubmitRequest(BaseModel):
    """주문 접수 요청"""
    order_id: str = Field(..., description="Platform에서 발급한 주문 ID")
    service_code: str = Field(..., description="서비스 코드 (nipt, carrier_screening, sgnipt)")
    sample_name: Optional[str] = Field(
        default=None,
        description="레거시 전용: 비어 있으면 파이프라인 폴더명은 order_id와 동일",
    )
    work_dir: Optional[str] = Field(default=None, description="작업 디렉토리명 (예: 2601)")
    
    # FASTQ 입력 (URL 또는 로컬 경로)
    fastq_r1_url: Optional[str] = Field(default=None, description="R1 FASTQ URL")
    fastq_r2_url: Optional[str] = Field(default=None, description="R2 FASTQ URL")
    fastq_r1_path: Optional[str] = Field(default=None, description="R1 FASTQ 로컬 경로")
    fastq_r2_path: Optional[str] = Field(default=None, description="R2 FASTQ 로컬 경로")
    
    # 서비스별 추가 파라미터
    params: Optional[Dict[str, Any]] = Field(
        default_factory=dict,
        description="서비스별 추가 파라미터 (예: backbone_bed, pon_tar 등)"
    )
    


class OrderStatusResponse(BaseModel):
    """주문 상태 응답"""
    order_id: str
    service_code: str
    status: OrderStatus
    progress: int = Field(default=0, ge=0, le=100)
    message: str = ""
    created_at: Optional[str] = None
    updated_at: Optional[str] = None


class OrderSubmitResponse(BaseModel):
    """주문 접수 응답"""
    status: str
    order_id: str
    service_code: str
    message: str
    queue_position: Optional[int] = None


class StartOrderRequest(BaseModel):
    """Force Run / Fresh Run 요청 본문 (선택)."""
    fresh: bool = Field(
        default=False,
        description="True면 파이프라인 캐시(Nextflow work/ 등)를 지우고 처음부터 재실행 (sgNIPT --fresh).",
    )
    use_ssd: bool = Field(
        default=False,
        description=(
            "Carrier / whole exome / health screening: run_analysis.sh 에 --use-ssd 및 "
            "(선택) --scratch-dir 를 붙여 SSD scratch 프로파일을 사용합니다."
        ),
    )
    scratch_dir: Optional[str] = Field(
        default=None,
        description="use_ssd 가 True일 때 run_analysis.sh --scratch-dir 에 넘길 호스트 경로 (예: /tmp/exome-scratch).",
    )


class OrderSaveResponse(BaseModel):
    """주문 저장만 (큐 미등록) 응답"""
    status: str
    order_id: str
    service_code: str
    message: str


class OrderUpdateRequest(BaseModel):
    """저장된·재실행 가능 주문 부분 수정 (PATCH)"""

    order_id: Optional[str] = Field(
        default=None,
        description="변경 시 새 주문 ID(URL 경로의 기존 ID는 삭제 후 이 키로 저장)",
    )
    sample_name: Optional[str] = None
    work_dir: Optional[str] = None
    fastq_r1_url: Optional[str] = None
    fastq_r2_url: Optional[str] = None
    fastq_r1_path: Optional[str] = None
    fastq_r2_path: Optional[str] = None
    params: Optional[Dict[str, Any]] = None


class OrderUpdateResponse(BaseModel):
    status: str
    order_id: str
    message: str


class UpdateFastqPathsRequest(BaseModel):
    """큐 대기 중인 주문의 로컬 FASTQ 경로 수정 (Portal 브라우저 선택)"""

    fastq_r1_path: Optional[str] = Field(default=None, description="R1 절대 경로; 빈 문자열이면 제거")
    fastq_r2_path: Optional[str] = Field(default=None, description="R2 절대 경로; 빈 문자열이면 제거")


class DarkGenesSectionReviewItem(BaseModel):
    """Index-aligned with ``dark_genes.detailed_sections`` (Portal Dark genes tab)."""

    approved: bool = Field(default=False, description="Include this block in PDF/report supplement")
    notes: str = Field(default="", description="Reviewer notes; shown in PDF when approved")
    risk: Optional[Literal["low", "high"]] = Field(
        default=None,
        description="PDF/portal title accent; omit to use pipeline default (low unless WARNING / warning-kind)",
    )


class DarkGenesReviewRequest(BaseModel):
    """Persist per-section approval, notes, and risk into ``result.json`` → ``dark_genes.section_reviews``."""

    section_reviews: List[DarkGenesSectionReviewItem] = Field(
        default_factory=list,
        description="One entry per detailed section, same order as detailed_sections",
    )


class PgxGeneReviewRow(BaseModel):
    """Per-gene review row aligned with ``pgx.gene_results[].gene`` (PharmCAT phenotype)."""

    gene: str = Field(..., min_length=1, max_length=64)
    reviewer_confirmed: bool = False
    reviewer_comment: str = Field(default="", max_length=4000)


class PgxCustomGeneReviewRow(BaseModel):
    """Per-variant review row for extended panel ``pgx.custom_gene_results[]`` (e.g. APOE rs429358 / rs7412)."""

    gene: str = Field(..., min_length=1, max_length=64)
    rsid: str = Field(..., min_length=1, max_length=32)
    reviewer_confirmed: bool = False
    reviewer_comment: str = Field(default="", max_length=4000)


class PgxReviewRequest(BaseModel):
    """Persist portal reviewer state into ``result.json`` → ``pgx.portal_review`` + per-gene fields on ``gene_results``."""

    reviewer_notes: str = Field(default="", max_length=16000, description="Internal notes; not in customer PDF by default")
    reviewed: bool = Field(default=False, description="Reviewer marked PGx section as reviewed")
    gene_reviews: List[PgxGeneReviewRow] = Field(
        default_factory=list,
        description="Updates reviewer_confirmed / reviewer_comment for matching gene symbols in pgx.gene_results",
    )
    custom_gene_reviews: List[PgxCustomGeneReviewRow] = Field(
        default_factory=list,
        description="Updates reviewer_confirmed / reviewer_comment for matching gene+rsid in pgx.custom_gene_results",
    )
    include_apoe_proactive_pdf: bool = Field(
        default=False,
        description="Include APOE ε2/ε3/ε4 proactive summary block in customer PDF (proactive health orders)",
    )


class WesPanelCustomSave(BaseModel):
    """Portal: create or replace a user-defined WES/exome panel (saved to ``wes_panels_custom_json``)."""

    id: str = Field(..., min_length=2, max_length=64, description="Stable slug (letters, digits, _ -)")
    label: str = Field(..., min_length=1, max_length=220)
    category: str = Field(default="other", max_length=64)
    description: Optional[str] = Field(default=None, max_length=4000)
    backbone_bed: Optional[str] = Field(default=None, description="Optional absolute or repo-relative backbone BED")
    disease_bed: Optional[str] = Field(
        default=None,
        description="Existing disease/panel BED path (absolute or repo-relative); used if genes not set",
    )
    genes: Optional[List[str]] = Field(
        default=None,
        description="HGNC symbols etc.; intervals taken from GENE_BED / settings.gene_bed",
    )
    genes_text: Optional[str] = Field(
        default=None,
        description="Alternative to genes: comma/newline-separated list",
    )
    gene_source_bed: Optional[str] = Field(
        default=None,
        description=(
            "BED to subset by gene symbols (e.g. Twist exome); column 4 = gene. "
            "Defaults to WES_PANEL_GENE_SOURCE_BED then GENE_BED."
        ),
    )
    skip_generated_bed: bool = Field(
        default=False,
        description=(
            "If true, only ``interpretation_genes`` are stored (post-analysis filtering). "
            "No generated panel BED; source interval BED is not required."
        ),
    )


# ─── Internal Job Model ───────────────────────────────────

class Job(BaseModel):
    """내부 작업 모델 (큐에서 관리)"""
    order_id: str
    service_code: str
    sample_name: str
    work_dir: str
    
    # 경로
    fastq_dir: Optional[str] = None
    analysis_dir: Optional[str] = None
    output_dir: Optional[str] = None
    log_dir: Optional[str] = None
    
    # 입력 소스
    fastq_r1_url: Optional[str] = None
    fastq_r2_url: Optional[str] = None
    fastq_r1_path: Optional[str] = None
    fastq_r2_path: Optional[str] = None
    
    # 서비스별 파라미터
    params: Dict[str, Any] = Field(default_factory=dict)
    
    # 상태
    status: OrderStatus = OrderStatus.RECEIVED
    progress: int = 0
    message: str = ""
    # 타임스탬프
    created_at: str = Field(default_factory=now_kst_iso)
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    updated_at: Optional[str] = None
    
    # 실행 정보
    pid: Optional[int] = None
    exit_code: Optional[int] = None
    duration: Optional[str] = None
    error_log: Optional[str] = None

    def update_status(self, status: OrderStatus, progress: int = None, message: str = None):
        """상태 업데이트"""
        self.status = status
        if progress is not None:
            self.progress = progress
        if message is not None:
            self.message = message
        self.updated_at = now_kst_iso()

    @model_validator(mode="after")
    def _normalize_carrier_legacy_container_paths(self) -> "Job":
        """
        gx-exome rename: SQLite job_json may still list /data/carrier_screening/... for paths
        and params (main_vcf, prior reuse hints). Map to /data/gx-exome* on load.
        """
        if self.service_code not in ("carrier_screening", "whole_exome", "health_screening"):
            return self
        from .config import normalize_legacy_carrier_container_path

        def norm(s: Optional[str]) -> Optional[str]:
            if not s or not isinstance(s, str):
                return s
            return normalize_legacy_carrier_container_path(s) or s

        updates: Dict[str, Any] = {}
        for fname in (
            "fastq_dir",
            "analysis_dir",
            "output_dir",
            "log_dir",
            "fastq_r1_path",
            "fastq_r2_path",
        ):
            v = getattr(self, fname)
            nv = norm(v)
            if nv != v:
                updates[fname] = nv

        p = dict(self.params or {})
        param_keys = (
            "main_vcf",
            "_prior_reuse_main_vcf_hint",
            "_prior_reuse_analysis_dir",
            "_prior_reuse_output_dir",
            "_prior_reuse_log_dir",
        )
        p_changed = False
        for k in param_keys:
            if k in p and isinstance(p[k], str):
                nk = norm(p[k])
                if nk != p[k]:
                    p[k] = nk
                    p_changed = True
        if p_changed:
            updates["params"] = p

        if updates:
            return self.model_copy(update=updates)
        return self


class NotificationResult(BaseModel):
    """알림 결과"""
    status: NotificationStatus
    message: str = ""
    response_code: Optional[int] = None


class OutputFile(BaseModel):
    """업로드할 출력 파일 정보"""
    file_path: str = Field(..., description="로컬 파일 경로")
    file_type: str = Field(..., description="파일 유형 (review_json, snapshot, vcf, pdf 등)")
    file_name: Optional[str] = Field(default=None, description="업로드 시 사용할 파일명")
    content_type: str = Field(default="application/octet-stream", description="MIME type")


class ReportGenerateRequest(BaseModel):
    """리포트 생성 요청 (리뷰어 확정 후 Portal에서 호출)"""
    confirmed_variants: List[Dict[str, Any]] = Field(
        ..., description="리뷰어가 확정한 변이 목록 (reviewer_confirmed, reviewer_classification, reviewer_comment 포함)"
    )
    reviewer_info: Dict[str, Any] = Field(
        ..., description="리뷰어 정보 (name, id, institution 등)"
    )
    patient_info: Optional[Dict[str, Any]] = Field(
        default=None,
        description="환자 정보 (name, dob, gender 등)"
    )
    partner_info: Optional[Dict[str, Any]] = Field(
        default=None,
        description="파트너 정보 (couple 검사 시, name, dob, gender 등)"
    )
    languages: Optional[List[str]] = Field(
        default=None,
        description="리포트 생성 언어 목록 (예: ['EN', 'CN', 'KO']). None이면 config 기본값 사용"
    )


class ReportGenerateResponse(BaseModel):
    """리포트 생성 응답"""
    status: str
    order_id: str
    service_code: str
    report_files: List[str] = Field(default_factory=list, description="생성된 리포트 파일 목록")
    uploaded_count: int = 0
    message: str = ""


class GeneKnowledgeSaveRequest(BaseModel):
    """Portal: persist edited gene_knowledge SQLite row (gene must appear in the order's variants)."""

    gene: str = Field(..., description="Gene symbol")
    function_summary: str = Field(default="", description="gene_data.function_summary")
    disease_association: str = Field(default="", description="gene_data.disease_association")
    disorder: str = Field(default="", description="gene_data.disorder")
    omim_number: str = Field(default="", description="gene_data.omim_number")
    inheritance: str = Field(default="", description="gene_data.inheritance")


class VariantKnowledgeSaveRequest(BaseModel):
    """Portal: per-variant row in ``variant_knowledge`` (stable key gene|HGVS; shared gene text stays in ``gene_data``)."""

    variant_key: str = Field(..., description="Same key as server make_variant_key(gene, hgvsc, hgvsp)")
    variant_notes: str = Field(default="", description="Variant-specific cache; repeats same key across orders")


class QueueSummary(BaseModel):
    """큐 상태 요약"""
    total_queued: int = 0
    total_running: int = 0
    total_completed: int = 0
    total_failed: int = 0
    jobs_by_service: Dict[str, int] = Field(
        default_factory=dict,
        description="서비스별 (queued + running) 합계 (하위 호환)",
    )
    stats_by_service: Dict[str, Dict[str, int]] = Field(
        default_factory=dict,
        description="서비스별 queued / running / completed / failed (QueueManager _stats)",
    )
    running_jobs: List[Dict[str, Any]] = Field(default_factory=list)
