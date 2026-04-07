"""
Service Daemon Configuration

환경변수 기반 설정 관리.
기존 nipt-daemon의 config.py를 일반화하여 여러 서비스를 지원합니다.
"""

import os
from typing import Optional, List
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field, field_validator, model_validator

# Renamed mount: compose uses /data/gx-exome; old .env may still reference /data/carrier_screening.
_LEGACY_CARRIER_DATA_ROOT = "/data/carrier_screening"
_GX_EXOME_DATA_ROOT = "/data/gx-exome"
_LEGACY_CARRIER_WORK_ROOT = "/data/carrier_screening_work"
_GX_EXOME_WORK_ROOT = "/data/gx-exome-work"


def normalize_legacy_carrier_container_path(value: Optional[str]) -> Optional[str]:
    """
    Map deprecated /data/carrier_screening* paths to gx-exome layout.

    Use for env-backed settings, persisted job.params paths (main_vcf, hints), and API payloads
    so SQLite rows from before the gx-exome rename still resolve on disk.
    """
    if not value or not isinstance(value, str):
        return value
    t = value.strip()
    if not t:
        return value
    if t == _LEGACY_CARRIER_DATA_ROOT or t.startswith(_LEGACY_CARRIER_DATA_ROOT + "/"):
        return _GX_EXOME_DATA_ROOT + t[len(_LEGACY_CARRIER_DATA_ROOT) :]
    if t == _LEGACY_CARRIER_WORK_ROOT or t.startswith(_LEGACY_CARRIER_WORK_ROOT + "/"):
        return _GX_EXOME_WORK_ROOT + t[len(_LEGACY_CARRIER_WORK_ROOT) :]
    return value


class Settings(BaseSettings):
    """서비스 데몬 설정"""

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # ─── Application ───────────────────────────────────────
    app_name: str = Field(default="service-daemon", description="Application name")
    app_env: str = Field(default="dev", description="Environment: dev, staging, prod")
    app_port: int = Field(default=8000, description="API server port")
    log_level: str = Field(default="INFO", description="Logging level")

    # ─── Inbound API protection (optional) ─────────────────
    # 비어 있으면 비활성. 설정 시 Authorization: Bearer <값> 또는 X-API-Key 헤더 필요 (portal/health/static 제외)
    api_key: Optional[str] = Field(
        default=None,
        description="Optional shared secret; enables Bearer / X-API-Key guard on API routes",
    )

    # ─── Platform API ──────────────────────────────────────
    platform_api_base: str = Field(
        default="https://api.genolyx.com",
        description="Platform API base URL"
    )
    platform_api_enabled: bool = Field(
        default=True,
        description=(
            "If false, skip Genolyx login and all platform calls (status, uploads, notify). "
            "Use for local/portal dev when credentials are not configured."
        ),
    )
    auth_url: Optional[str] = Field(
        default=None,
        description="Authentication endpoint URL"
    )
    api_username: Optional[str] = Field(default=None, description="API login username")
    api_password: Optional[str] = Field(default=None, description="API login password")

    # ─── Telegram (optional job lifecycle alerts) ──────────
    telegram_bot_token: Optional[str] = Field(
        default=None,
        description="Telegram Bot API token; unset disables sendMessage notifications",
    )
    telegram_chat_ids: Optional[str] = Field(
        default=None,
        description="Comma-separated chat IDs (e.g. personal DM or group)",
    )
    telegram_notify_enabled: bool = Field(
        default=True,
        description="If false, do not send Telegram messages",
    )

    # ─── Queue / Worker ────────────────────────────────────
    max_concurrent_jobs: int = Field(
        default=2,
        description="Maximum number of concurrent pipeline jobs"
    )
    queue_poll_interval: int = Field(
        default=30,
        description="Seconds between queue polling cycles"
    )

    # ─── Base Directories ──────────────────────────────────
    base_dir: str = Field(
        default="/data",
        description="Root data directory"
    )
    orders_db_path: Optional[str] = Field(
        default=None,
        description=(
            "SQLite file for order/job persistence and JSON snapshots "
            "(result, review, report). Default: {BASE_DIR}/service-daemon/orders.db"
        ),
    )
    fastq_base_dir: str = Field(
        default="/data/fastq",
        description="Default FASTQ directory (nipt 등, 포털 browse 미지정 시)"
    )
    sgnipt_fastq_dir: str = Field(
        default="/home/ken/sgNIPT/fastq",
        description="sgNIPT: 포털 FASTQ browse 루트 (비우면 sgnipt_job_root/fastq)"
    )
    sgnipt_data_dir: str = Field(
        default="/home/ken/sgNIPT/data",
        description="sgNIPT: BAM samplesheet CSV browse 루트 (비우면 sgnipt_job_root/data)"
    )
    sgnipt_layout_root: str = Field(
        default="/home/ken/sgNIPT",
        description=(
            "sgNIPT 논리 루트(참조). 실제 주문 경로는 sgnipt_work_root 가 있으면 그쪽을 사용."
        ),
    )
    sgnipt_work_root: Optional[str] = Field(
        default=None,
        description=(
            "sgNIPT 쓰기 가능 루트: fastq|analysis|output|log/<work_dir>/<order_id>/ "
            "미설정 시 sgnipt_layout_root 사용(해당 경로가 데몬 UID에 쓰기 가능해야 함)."
        ),
    )
    sgnipt_docker_image: str = Field(
        default="sgnipt",
        description="sgNIPT Docker 이미지명 (docker run)",
    )
    sgnipt_run_script_path: str = Field(
        default="/home/ken/sgNIPT/src/run_sgnipt.sh",
        description="호스트에서 bash 로 실행할 run_sgnipt.sh 절대 경로 (레포 루트의 src/ 아래)",
    )
    sgnipt_src_root: Optional[str] = Field(
        default=None,
        description=(
            "sgNIPT 소스 클론 루트(선택). Docker 를 풀린 레포(/home/ken/sgnipt) 와 데이터 트리(/home/ken/sgNIPT) 가 "
            "다를 때 run_sgnipt.sh 탐색용. 비우면 SGNIPT_RUN_SCRIPT_PATH · WORK/LAYOUT 만 시도."
        ),
    )
    sgnipt_container_mount_root: str = Field(
        default="/Work/SgNIPT",
        description="docker run -v <job_root>:<이 경로> (sgnipt 이미지의 /Work/SgNIPT)",
    )
    sgnipt_use_docker: bool = Field(
        default=True,
        description="True면 docker run, False면 호스트에서 bash 직접 실행",
    )
    sgnipt_docker_extra_args: str = Field(
        default="",
        description="docker run --rm 뒤 이미지 앞에 넣는 추가 인자(공백 구분)",
    )
    carrier_screening_fastq_dir: str = Field(
        default="/home/ken/gx-exome/fastq",
        description="Carrier screening: 포털 FASTQ browse 및 경로 검증용 루트"
    )
    analysis_base_dir: str = Field(
        default="/data/analysis",
        description="Base directory for analysis working directories"
    )
    output_base_dir: str = Field(
        default="/data/output",
        description="Base directory for final output files"
    )
    log_base_dir: str = Field(
        default="/data/log",
        description="Base directory for pipeline logs"
    )

    # ─── Service-specific Directories ──────────────────────
    # 각 서비스별 파이프라인 경로 (환경변수로 오버라이드 가능)
    carrier_screening_pipeline_dir: str = Field(
        default="/opt/pipelines/carrier-screening",
        description="Carrier Screening Nextflow pipeline directory"
    )
    carrier_screening_main_nf: Optional[str] = Field(
        default=None,
        description=(
            "Absolute path to main.nf if not at pipeline_dir/bin/main.nf "
            "(e.g. /opt/pipelines/carrier-screening/main.nf). "
            "Can also point at dark_gene_pipeline/main.nf (see dark_gene_pipeline_dir)."
        ),
    )
    dark_gene_pipeline_dir: Optional[str] = Field(
        default=None,
        description=(
            "Optional path to the unified dark-gene Nextflow repo (main.nf: FASTQ→align→dark-genes tracks). "
            "When set and main.nf exists, the daemon prefers this entry over CARRIER_SCREENING_PIPELINE_DIR "
            "unless CARRIER_SCREENING_MAIN_NF is set. "
            "Use CARRIER_SCREENING_FORCE_DIRECT_NEXTFLOW=true (or remove run_analysis.sh) so jobs use "
            "nextflow run directly; otherwise auto-discovered run_analysis.sh takes precedence."
        ),
    )
    dark_gene_skip_cnv: bool = Field(
        default=True,
        description="When using dark_gene unified main.nf, pass --skip_cnv (gCNV is heavy; set false to enable).",
    )
    carrier_screening_force_direct_nextflow: bool = Field(
        default=False,
        description=(
            "When true, never use run_analysis.sh (even if present under layout or pipeline dir); "
            "run `nextflow run` with DARK_GENE_PIPELINE_DIR / CARRIER_SCREENING_MAIN_NF resolution instead."
        ),
    )
    carrier_screening_nextflow_profile: Optional[str] = Field(
        default=None,
        description=(
            "Optional Nextflow profile (e.g. docker, singularity) appended as -profile <name> "
            "for carrier direct Nextflow runs (not run_analysis.sh)."
        ),
    )
    carrier_screening_run_script: Optional[str] = Field(
        default=None,
        description=(
            "Optional: host path to run_analysis.sh (Dark Gene wrapper). "
            "If set, carrier jobs run this script instead of direct nextflow; requires Docker socket."
        ),
    )
    carrier_screening_script_data_dir: Optional[str] = Field(
        default=None,
        description="-d for run_analysis: project root with fastq/, data/refs, data/bed (default: BASE_DIR/carrier-screening)",
    )
    carrier_screening_script_ref_dir: Optional[str] = Field(
        default=None,
        description="-r for run_analysis (default: <data_dir>/data/refs)",
    )
    carrier_screening_script_extra_args: str = Field(
        default="",
        description=(
            "Additional args for run_analysis.sh (shlex-split, appended last). "
            "--work-dir/--sample/--data-dir/--panel/--aligner/--variant-caller/--no-skip-vep/--skip-cnv "
            "are already built explicitly; add only truly extra flags here (e.g. --shared-ref-dir /path)."
        ),
    )
    carrier_screening_artifact_base: Optional[str] = Field(
        default=None,
        description=(
            "Writable root for carrier fastq staging and analysis/output/log (subdirs fastq|analysis|output|log). "
            "When unset, those paths use carrier_screening_layout_base. "
            "Set when the layout tree is not writable by the daemon UID (e.g. /data/gx-exome-work)."
        ),
    )
    carrier_screening_report_output_root: Optional[str] = Field(
        default=None,
        description=(
            "If set, Portal Generate Report writes report.json and PDFs under "
            "<root>/output/<work_dir>/<sample_name> and reads QC from "
            "<root>/analysis/<work_dir>/<sample_name> (same layout as pipeline output). "
            "Unset uses carrier artifact/output paths as today. Example: /home/sam/Carrier_result"
        ),
    )
    carrier_screening_report_template_dir: Optional[str] = Field(
        default=None,
        description=(
            "Jinja2 HTML templates for carrier PDF (carrier_EN.html, carrier_couples_EN.html, …). "
            "Highest priority. Relative paths are resolved from the repo root (not process cwd). "
            "If unset: <CARRIER_SCREENING_REPORT_OUTPUT_ROOT>/report_templates (or legacy "
            "…/carrier_report), then packaged data/report_templates, then REPORT_TEMPLATE_DIR, "
            "then pipeline data/templates."
        ),
    )
    carrier_capture_panel_bed_dir: Optional[str] = Field(
        default=None,
        description=(
            "Capture panel BED root for run_analysis.sh --backbone-bed. "
            "Layout: {dir}/{capture_panel_id}/targets.bed(.gz/.gz.tbi). "
            "Must be a HOST path accessible inside the gx-exome inner docker container "
            "(i.e. under CARRIER_SCREENING_SCRIPT_DATA_DIR/data/bed, or any path "
            "explicitly mounted in run_analysis.sh docker run). "
            "Default: {CARRIER_SCREENING_SCRIPT_DATA_DIR}/data/bed."
        ),
    )
    carrier_default_backbone_bed: Optional[str] = Field(
        default=None,
        description="Default backbone BED if job params omit backbone_bed (portal override still wins)",
    )
    carrier_default_disease_bed: Optional[str] = Field(
        default=None,
        description="Default disease/panel BED if job params omit disease_bed",
    )
    wes_panels_json: Optional[str] = Field(
        default=None,
        description=(
            "Optional path to WES/exome panel catalog JSON (default: data/wes_panels.json under repo root). "
            "Relative paths resolve from repo root."
        ),
    )
    wes_panels_custom_json: Optional[str] = Field(
        default=None,
        description=(
            "Portal-saved custom panels (default: data/wes_panels_custom.json). "
            "Merged with wes_panels_json; custom entries override same id."
        ),
    )
    wes_panels_generated_dir: Optional[str] = Field(
        default=None,
        description=(
            "Directory for BED files built from gene lists in the panel builder "
            "(default: data/bed/wes_panels_generated under repo root)."
        ),
    )
    wes_panel_gene_source_bed: Optional[str] = Field(
        default=None,
        description=(
            "Default BED used when building a panel from a gene list (Twist/exome capture, etc.); "
            "column 4 = gene symbol. Falls back to GENE_BED if unset."
        ),
    )
    nipt_pipeline_dir: str = Field(
        default="/opt/pipelines/nipt",
        description="NIPT pipeline directory"
    )
    sgnipt_pipeline_dir: str = Field(
        default="/opt/pipelines/sgnipt",
        description="sgNIPT pipeline directory"
    )

    # ─── Nextflow ──────────────────────────────────────────
    nextflow_executable: str = Field(
        default="nextflow",
        description="Path to Nextflow executable"
    )
    nextflow_config: Optional[str] = Field(
        default=None,
        description="Path to global Nextflow config file"
    )

    # ─── Reference Data ────────────────────────────────────
    ref_fasta: Optional[str] = Field(default=None, description="Reference genome FASTA")
    ref_fai: Optional[str] = Field(default=None, description="Reference genome FAI index")
    ref_dict: Optional[str] = Field(default=None, description="Reference genome dict")
    ref_bwa_indices: Optional[str] = Field(default=None, description="BWA index directory")

    # ─── Annotation Resources ──────────────────────────────
    clinvar_vcf: Optional[str] = Field(default=None, description="ClinVar VCF path")
    gnomad_vcf: Optional[str] = Field(default=None, description="gnomAD VCF path")
    gnomad_dir: Optional[str] = Field(default=None, description="gnomAD VCF directory (per-chromosome)")
    gnomad_genomes_glob: str = Field(
        default="gnomad.genomes.v*.sites*.bgz",
        description="Glob pattern for gnomAD genomes VCF files"
    )
    gnomad_exomes_glob: str = Field(
        default="gnomad.exomes.v*.sites*.bgz",
        description="Glob pattern for gnomAD exomes VCF files"
    )
    dbsnp_vcf: Optional[str] = Field(default=None, description="dbSNP VCF path")
    snpeff_jar: Optional[str] = Field(default=None, description="snpEff JAR path")
    snpeff_db: str = Field(default="GRCh38.86", description="snpEff genome DB name (snpEff download -v)")
    snpeff_data_dir: Optional[str] = Field(
        default=None,
        description="snpEff -dataDir (persistent genomes; e.g. /data/reference/snpeff on host)",
    )

    # ─── ClinGen ───────────────────────────────────────────
    clingen_tsv: Optional[str] = Field(
        default=None,
        description="ClinGen gene curation TSV path"
    )

    # ─── MANE RefSeq ───────────────────────────────────────
    mane_gff: Optional[str] = Field(
        default=None,
        description="MANE RefSeq GFF file path"
    )

    # ─── Gene Coordinates ──────────────────────────────────
    gene_bed: Optional[str] = Field(
        default=None,
        description="Gene coordinates BED file (genes_hg38.bed)"
    )

    # ─── HPO (Human Phenotype Ontology) ────────────────────
    hpo_gene_file: Optional[str] = Field(
        default=None,
        description="HPO genes_to_phenotype.txt file path"
    )

    # ─── Curated Variant Database ──────────────────────────
    curated_variants_db: Optional[str] = Field(
        default=None,
        description="SQLite database path for curated variant classifications"
    )

    # ─── Disease Database ──────────────────────────────────
    disease_db_dir: Optional[str] = Field(
        default=None,
        description="Directory containing disease-gene mapping databases"
    )
    disease_gene_json: Optional[str] = Field(
        default=None,
        description="JSON file mapping diseases to genes and inheritance patterns"
    )
    gene_knowledge_db: Optional[str] = Field(
        default=None,
        description=(
            "SQLite path (gene_data table). Used only when generating a report for "
            "reviewer-confirmed variants — not during bulk VCF annotation. "
            "If unset, may default to genetic_reporter_lookup's genetic_knowledge.db (see gene_knowledge_db_fallback_path)."
        ),
    )
    gene_knowledge_db_fallback_path: Optional[str] = Field(
        default="/home/sam/genetic_reporter_lookup/genetic_knowledge.db",
        description=(
            "When gene_knowledge_db is empty: if the parent directory of this path exists, use it as "
            "gene_knowledge_db (same file and gene_data schema as genetic_reporter_lookup app2.py)."
        ),
    )
    gene_knowledge_enrich_on_report: bool = Field(
        default=True,
        description=(
            "If True and gene_knowledge_db is set, merge disorder/inheritance from that DB "
            "(and optionally Gemini) into confirmed variants during report generation."
        ),
    )
    gene_knowledge_gemini_on_report: bool = Field(
        default=True,
        description=(
            "If True, gene_knowledge_db cache miss during report enrichment may call Gemini "
            "(requires gemini_api_key). Set False to use SQLite cache only."
        ),
    )
    gene_knowledge_gemini_model: str = Field(
        default="gemini-2.5-flash",
        description="Gemini model for gene_knowledge_db Search+JSON extraction.",
    )

    # ─── HGMD (Human Gene Mutation Database) ───────────────
    hgmd_vcf: Optional[str] = Field(
        default=None,
        description="HGMD Professional VCF path (licensed)"
    )

    # ─── Known Gene Lists ──────────────────────────────────
    known_carrier_genes: Optional[str] = Field(
        default=None,
        description="File containing known carrier screening gene list"
    )
    acmg_sf_genes: Optional[str] = Field(
        default=None,
        description="ACMG Secondary Findings gene list (v3.2+)"
    )

    # ─── Report Templates ──────────────────────────────────
    report_template_dir: Optional[str] = Field(
        default=None,
        description=(
            "Legacy: Jinja2 HTML report templates. For PDFs, packaged data/report_templates "
            "is preferred when set; use carrier_screening_report_template_dir to override. "
            "Relative paths are resolved from the repo root."
        ),
    )
    report_languages: str = Field(
        default="EN,CN",
        description="Comma-separated list of report languages"
    )
    report_logo_path: Optional[str] = Field(
        default=None,
        description="Path to company logo for reports"
    )

    # ─── AI Classification ─────────────────────────────────
    acmg_ai_enabled: bool = Field(
        default=False,
        description="Enable AI-assisted ACMG classification"
    )
    acmg_ai_provider: str = Field(
        default="gemini",
        description="AI provider: gemini | openai"
    )
    acmg_ai_model: str = Field(
        default="gemini-2.5-flash",
        description="AI model name (e.g. gemini-2.5-flash or gpt-4.1-mini)"
    )
    gemini_api_key: str = Field(
        default="",
        description="Google Gemini API key (env GEMINI_API_KEY; if unset, may load from gemini_api_key_env_file)",
    )
    gemini_api_key_env_file: Optional[str] = Field(
        default="/home/sam/genetic_reporter_lookup/.env",
        description=(
            "If GEMINI_API_KEY is empty, load it from this .env (same convention as genetic_reporter_lookup). "
            "Override with env GEMINI_API_KEY_ENV_FILE or set to empty to disable."
        ),
    )
    acmg_ai_api_key: str = Field(
        default="",
        description="OpenAI API key (used when acmg_ai_provider=openai)"
    )

    # ─── Literature Search (PubMed) ────────────────────────
    literature_enabled: bool = Field(
        default=True,
        description="Enable PubMed literature search for variants (cached permanently)"
    )
    literature_db_path: Optional[str] = Field(
        default=None,
        description="SQLite path for literature cache (default: {base_dir}/service-daemon/literature.db)"
    )
    literature_max_results: int = Field(
        default=10,
        description="Maximum PubMed articles to fetch per variant query"
    )
    ncbi_email: str = Field(
        default="service-daemon@example.com",
        description="NCBI E-utilities contact email (required by NCBI policy)"
    )
    ncbi_api_key: str = Field(
        default="",
        description="NCBI API key (optional; increases rate limit from 3 to 10 req/s)"
    )
    ncbi_tool: str = Field(
        default="service-daemon",
        description="NCBI E-utilities tool name"
    )

    # ─── Enabled Services ──────────────────────────────────
    enabled_services: str = Field(
        default="carrier_screening,whole_exome,health_screening,nipt,sgnipt",
        description=(
            "Comma-separated list of enabled service codes. "
            "Must include any code the portal exposes (e.g. whole_exome, health_screening) "
            "or save/submit will return Unknown service_code."
        ),
    )

    @field_validator("gemini_api_key", "acmg_ai_api_key", mode="before")
    @classmethod
    def _strip_api_key_fields(cls, v: object) -> object:
        if isinstance(v, str):
            return v.strip()
        return v

    @property
    def resolved_orders_db_path(self) -> str:
        """Effective SQLite path for orders (always under base_dir unless overridden)."""
        if self.orders_db_path and str(self.orders_db_path).strip():
            return os.path.abspath(str(self.orders_db_path).strip())
        return os.path.join(self.base_dir, "service-daemon", "orders.db")

    @property
    def resolved_literature_db_path(self) -> str:
        """Effective SQLite path for literature cache."""
        if self.literature_db_path and str(self.literature_db_path).strip():
            return os.path.abspath(str(self.literature_db_path).strip())
        return os.path.join(self.base_dir, "service-daemon", "literature.db")

    @property
    def sgnipt_job_root(self) -> str:
        """주문별 fastq|analysis|output|log 의 상위 루트 (쓰기 가능 경로 권장)."""
        wr = (self.sgnipt_work_root or "").strip()
        if wr:
            return os.path.abspath(wr)
        lr = (self.sgnipt_layout_root or "").strip() or "/home/ken/sgNIPT"
        return os.path.abspath(lr)

    @property
    def sgnipt_fastq_root(self) -> str:
        """FASTQ browse/검증 루트 (명시 fastq_dir 또는 job_root/fastq)."""
        fq = (self.sgnipt_fastq_dir or "").strip()
        if fq:
            return os.path.abspath(fq)
        return os.path.abspath(os.path.join(self.sgnipt_job_root, "fastq"))

    @property
    def carrier_screening_layout_base(self) -> str:
        """
        Carrier FASTQ 트리 상위 (fastq/ 가 그 아래).
        CARRIER_SCREENING_FASTQ_DIR 가 .../fastq 로 끝나면 그 부모를 씀 (예: /data/gx-exome).
        """
        fq = (self.carrier_screening_fastq_dir or "").strip().rstrip("/")
        if fq.lower().endswith("/fastq"):
            return os.path.abspath(os.path.dirname(fq))
        return os.path.join(self.base_dir, "gx-exome")

    @property
    def carrier_screening_work_root(self) -> str:
        """analysis / output / log 를 둘 쓰기 가능 루트 (미설정 시 layout_base 와 동일)."""
        ab = (self.carrier_screening_artifact_base or "").strip()
        if ab:
            return os.path.abspath(ab)
        return self.carrier_screening_layout_base

    @property
    def enabled_service_list(self) -> List[str]:
        """활성화된 서비스 코드 목록"""
        return [s.strip() for s in self.enabled_services.split(",") if s.strip()]

    @property
    def report_language_list(self) -> List[str]:
        """리포트 언어 목록"""
        return [s.strip() for s in self.report_languages.split(",") if s.strip()]

    def get_carrier_screening_data_dir(self) -> str:
        """Carrier Screening 파이프라인의 data 디렉토리 경로"""
        return os.path.join(self.carrier_screening_pipeline_dir, "data")

    def get_carrier_screening_bed_dir(self) -> str:
        """Carrier Screening BED 파일 디렉토리 경로"""
        return os.path.join(self.carrier_screening_pipeline_dir, "data", "bed")

    @model_validator(mode="after")
    def _fill_gemini_api_key_from_genetic_reporter_lookup(self) -> "Settings":
        """Match genetic_reporter_lookup: GEMINI_API_KEY from that project's .env when not set here."""
        if (self.gemini_api_key or "").strip():
            return self
        path = (self.gemini_api_key_env_file or "").strip()
        if not path or not os.path.isfile(path):
            return self
        try:
            from dotenv import dotenv_values

            vals = dotenv_values(path)
            k = (vals.get("GEMINI_API_KEY") or "").strip()
            if k:
                return self.model_copy(update={"gemini_api_key": k})
        except Exception:
            pass
        return self

    @model_validator(mode="after")
    def _fill_gene_knowledge_db_from_genetic_reporter_lookup(self) -> "Settings":
        """Use the same SQLite as genetic_reporter_lookup (genetic_knowledge.db) when not set."""
        if (self.gene_knowledge_db or "").strip():
            return self
        path = (self.gene_knowledge_db_fallback_path or "").strip()
        if not path:
            return self
        abspath = os.path.abspath(path)
        parent = os.path.dirname(abspath)
        if not os.path.isdir(parent):
            return self
        return self.model_copy(update={"gene_knowledge_db": abspath})

    @field_validator(
        "carrier_screening_fastq_dir",
        "carrier_screening_run_script",
        "carrier_screening_script_data_dir",
        "carrier_screening_script_ref_dir",
        "carrier_screening_artifact_base",
        "carrier_screening_report_output_root",
        mode="before",
    )
    @classmethod
    def _normalize_legacy_carrier_paths(cls, v: object) -> object:
        """Map /data/carrier_screening* → /data/gx-exome* when env still uses old mount paths."""
        if v is None:
            return v
        if not isinstance(v, str):
            return v
        return normalize_legacy_carrier_container_path(v)


# 전역 설정 인스턴스
settings = Settings()
