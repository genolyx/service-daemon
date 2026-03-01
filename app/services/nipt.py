"""
NIPT Service Plugin

기존 nipt-daemon의 파이프라인 실행 로직을 플러그인으로 이전합니다.
"""

import os
import glob
import logging
from typing import List, Dict, Any

from .base import ServicePlugin
from app.config import settings
from app.models import Job, OutputFile

logger = logging.getLogger(__name__)


class NIPTPlugin(ServicePlugin):
    """NIPT 서비스 플러그인"""

    @property
    def service_code(self) -> str:
        return "nipt"

    @property
    def display_name(self) -> str:
        return "NIPT (Non-Invasive Prenatal Testing)"

    async def prepare_inputs(self, job: Job) -> bool:
        """FASTQ 파일 준비"""
        try:
            fastq_dir = job.fastq_dir
            os.makedirs(fastq_dir, exist_ok=True)

            # URL이 제공된 경우 다운로드
            if job.fastq_r1_url and job.fastq_r2_url:
                import subprocess
                r1_path = os.path.join(fastq_dir, os.path.basename(job.fastq_r1_url))
                r2_path = os.path.join(fastq_dir, os.path.basename(job.fastq_r2_url))

                for url, path in [(job.fastq_r1_url, r1_path), (job.fastq_r2_url, r2_path)]:
                    if not os.path.exists(path):
                        logger.info(f"[nipt] Downloading {url} -> {path}")
                        subprocess.run(
                            ["wget", "-q", "-O", path, url],
                            check=True, timeout=3600
                        )

            # 로컬 경로가 제공된 경우 심볼릭 링크
            elif job.fastq_r1_path and job.fastq_r2_path:
                for src in [job.fastq_r1_path, job.fastq_r2_path]:
                    dst = os.path.join(fastq_dir, os.path.basename(src))
                    if not os.path.exists(dst):
                        os.symlink(src, dst)

            # FASTQ 파일 존재 확인
            fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq.gz")) + \
                          glob.glob(os.path.join(fastq_dir, "*.fq.gz"))
            if not fastq_files:
                logger.error(f"[nipt] No FASTQ files found in {fastq_dir}")
                return False

            logger.info(f"[nipt] Input ready: {len(fastq_files)} FASTQ files in {fastq_dir}")
            return True

        except Exception as e:
            logger.error(f"[nipt] Failed to prepare inputs: {e}")
            return False

    async def get_pipeline_command(self, job: Job) -> str:
        """NIPT 파이프라인 실행 명령어"""
        pipeline_dir = settings.nipt_pipeline_dir
        analysis_dir = job.analysis_dir
        os.makedirs(analysis_dir, exist_ok=True)

        # 기존 nipt-daemon의 runner.py에서 사용하던 Docker 기반 명령어 구조
        cmd = (
            f"cd {pipeline_dir} && "
            f"bash run_nipt.sh "
            f"--fastq_dir {job.fastq_dir} "
            f"--output_dir {analysis_dir} "
            f"--sample_name {job.sample_name}"
        )

        # 추가 파라미터
        for key, value in job.params.items():
            cmd += f" --{key} {value}"

        return cmd

    async def check_completion(self, job: Job) -> bool:
        """NIPT 파이프라인 완료 확인"""
        analysis_dir = job.analysis_dir
        
        # 결과 파일 확인 (NIPT 특화)
        result_files = glob.glob(os.path.join(analysis_dir, "**", "*result*.json"), recursive=True)
        if result_files:
            logger.info(f"[nipt] Found result files: {result_files}")
            return True

        logger.warning(f"[nipt] No result files found in {analysis_dir}")
        return False

    async def process_results(self, job: Job) -> bool:
        """NIPT 결과 후처리 (PDF 생성 등)"""
        try:
            output_dir = job.output_dir
            os.makedirs(output_dir, exist_ok=True)

            # 기존 nipt-daemon의 리포트 생성 로직 호출
            # (실제 구현은 기존 nipt-daemon의 report 모듈에 의존)
            logger.info(f"[nipt] Processing results for {job.order_id}")
            
            # 결과 파일을 output_dir로 복사
            import shutil
            analysis_dir = job.analysis_dir
            for pattern in ["*result*.json", "*report*.pdf", "*.png"]:
                for f in glob.glob(os.path.join(analysis_dir, "**", pattern), recursive=True):
                    dst = os.path.join(output_dir, os.path.basename(f))
                    shutil.copy2(f, dst)
                    logger.info(f"[nipt] Copied {f} -> {dst}")

            return True

        except Exception as e:
            logger.error(f"[nipt] Result processing failed: {e}")
            return False

    async def get_output_files(self, job: Job) -> List[OutputFile]:
        """업로드할 결과 파일 목록"""
        output_dir = job.output_dir
        files = []

        # JSON 결과
        for f in glob.glob(os.path.join(output_dir, "*.json")):
            files.append(OutputFile(
                file_path=f,
                file_type="result_json",
                content_type="application/json"
            ))

        # PDF 리포트
        for f in glob.glob(os.path.join(output_dir, "*.pdf")):
            files.append(OutputFile(
                file_path=f,
                file_type="pdf_report",
                content_type="application/pdf"
            ))

        # 이미지
        for f in glob.glob(os.path.join(output_dir, "*.png")):
            files.append(OutputFile(
                file_path=f,
                file_type="image",
                content_type="image/png"
            ))

        logger.info(f"[nipt] Output files: {len(files)} files for upload")
        return files

    def get_progress_stages(self) -> Dict[str, int]:
        """NIPT 파이프라인 진행률 단계"""
        return {
            "Alignment": 20,
            "Deduplication": 30,
            "GC_Correction": 50,
            "Z_Score": 70,
            "Report_Generation": 90,
            "Completed": 100
        }


def create_plugin() -> ServicePlugin:
    """플러그인 팩토리 함수"""
    return NIPTPlugin()
