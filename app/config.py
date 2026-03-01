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
    dbsnp_vcf: Optional[str] = Field(default=None, description="dbSNP VCF path")
    snpeff_jar: Optional[str] = Field(default=None, description="snpEff JAR path")
    snpeff_db: str = Field(default="GRCh38.86", description="snpEff database name")

    # ─── Enabled Services ──────────────────────────────────
    enabled_services: str = Field(
        default="carrier_screening,nipt,sgnipt",
        description="Comma-separated list of enabled service codes"
    )

    @property
    def enabled_service_list(self) -> List[str]:
        """활성화된 서비스 코드 목록"""
        return [s.strip() for s in self.enabled_services.split(",") if s.strip()]

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"
        extra = "ignore"


# 전역 설정 인스턴스
settings = Settings()
