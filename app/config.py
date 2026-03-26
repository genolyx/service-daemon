"""
Service Daemon Configuration

환경변수 기반 설정 관리.
기존 nipt-daemon의 config.py를 일반화하여 여러 서비스를 지원합니다.
"""

import os
from typing import Optional, List
from pydantic_settings import BaseSettings
from pydantic import Field


class Settings(BaseSettings):
    """서비스 데몬 설정"""

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
        description="sgNIPT: 포털 FASTQ browse 및 경로 검증용 루트"
    )
    carrier_screening_fastq_dir: str = Field(
        default="/home/ken/carrier_screening/fastq",
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
            "(e.g. /opt/pipelines/carrier-screening/main.nf)"
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
        default="--skip-cnv",
        description="Extra args for run_analysis.sh (e.g. empty string to allow CNV)",
    )
    carrier_screening_artifact_base: Optional[str] = Field(
        default=None,
        description=(
            "Writable root for carrier analysis/output/log (subdirs analysis|output|log). "
            "FASTQ paths stay under carrier_screening_layout_base. "
            "Set when project/FASTQ tree is not writable by the container UID (e.g. use /data/carrier_screening_work)."
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
        description="Directory containing Jinja2 HTML report templates"
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
        default="gemini-2.0-flash",
        description="AI model name (e.g. gemini-2.0-flash or gpt-4.1-mini)"
    )
    gemini_api_key: str = Field(
        default="",
        description="Google Gemini API key"
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
        default="carrier_screening,nipt,sgnipt",
        description="Comma-separated list of enabled service codes"
    )

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
    def carrier_screening_layout_base(self) -> str:
        """
        Carrier FASTQ 트리 상위 (fastq/ 가 그 아래).
        CARRIER_SCREENING_FASTQ_DIR 가 .../fastq 로 끝나면 그 부모를 씀 (예: /data/carrier_screening).
        """
        fq = (self.carrier_screening_fastq_dir or "").strip().rstrip("/")
        if fq.lower().endswith("/fastq"):
            return os.path.abspath(os.path.dirname(fq))
        return os.path.join(self.base_dir, "carrier-screening")

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

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"
        extra = "ignore"


# 전역 설정 인스턴스
settings = Settings()
