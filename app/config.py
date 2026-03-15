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

    # ─── Platform API ──────────────────────────────────────
    platform_api_base: str = Field(
        default="https://api.genolyx.com",
        description="Platform API base URL"
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
    fastq_base_dir: str = Field(
        default="/data/fastq",
        description="Base directory for FASTQ files"
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
    snpeff_db: str = Field(default="GRCh38.86", description="snpEff database name")

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
    acmg_ai_model: str = Field(
        default="gpt-4.1-mini",
        description="AI model for ACMG classification"
    )

    # ─── Enabled Services ──────────────────────────────────
    enabled_services: str = Field(
        default="carrier_screening,nipt,sgnipt",
        description="Comma-separated list of enabled service codes"
    )

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
