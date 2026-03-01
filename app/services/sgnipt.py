"""
sgNIPT Service Plugin (Stub)

sgNIPT 서비스의 플러그인 스텁입니다.
실제 파이프라인 로직은 추후 구현합니다.
"""

import os
import glob
import logging
from typing import List, Dict, Any

from .base import ServicePlugin
from app.config import settings
from app.models import Job, OutputFile

logger = logging.getLogger(__name__)


class SgNIPTPlugin(ServicePlugin):
    """sgNIPT 서비스 플러그인 (스텁)"""

    @property
    def service_code(self) -> str:
        return "sgnipt"

    @property
    def display_name(self) -> str:
        return "sgNIPT (Single-Gene NIPT)"

    async def prepare_inputs(self, job: Job) -> bool:
        """FASTQ 파일 준비"""
        try:
            os.makedirs(job.fastq_dir, exist_ok=True)
            os.makedirs(job.analysis_dir, exist_ok=True)
            os.makedirs(job.output_dir, exist_ok=True)

            if job.fastq_r1_path and job.fastq_r2_path:
                for src in [job.fastq_r1_path, job.fastq_r2_path]:
                    dst = os.path.join(job.fastq_dir, os.path.basename(src))
                    if not os.path.exists(dst):
                        os.symlink(os.path.abspath(src), dst)

            logger.info(f"[sgnipt] Input prepared for {job.order_id}")
            return True
        except Exception as e:
            logger.error(f"[sgnipt] Failed to prepare inputs: {e}")
            return False

    async def get_pipeline_command(self, job: Job) -> str:
        """sgNIPT 파이프라인 실행 명령어"""
        pipeline_dir = settings.sgnipt_pipeline_dir
        cmd = (
            f"cd {pipeline_dir} && "
            f"bash run_sgnipt.sh "
            f"--fastq_dir {job.fastq_dir} "
            f"--output_dir {job.analysis_dir} "
            f"--sample_name {job.sample_name}"
        )
        return cmd

    async def check_completion(self, job: Job) -> bool:
        """파이프라인 완료 확인"""
        result_files = glob.glob(
            os.path.join(job.analysis_dir, "**", "*result*"), recursive=True
        )
        return len(result_files) > 0

    async def process_results(self, job: Job) -> bool:
        """결과 후처리"""
        logger.info(f"[sgnipt] Processing results for {job.order_id}")
        # TODO: sgNIPT 결과 후처리 구현
        return True

    async def get_output_files(self, job: Job) -> List[OutputFile]:
        """출력 파일 목록"""
        files = []
        for f in glob.glob(os.path.join(job.output_dir, "**", "*"), recursive=True):
            if os.path.isfile(f):
                files.append(OutputFile(
                    file_path=f,
                    file_type="result",
                    content_type="application/octet-stream"
                ))
        return files

    def get_progress_stages(self) -> Dict[str, int]:
        return {
            "Alignment": 20,
            "Variant_Calling": 50,
            "Analysis": 80,
            "Report": 95
        }


def create_plugin() -> ServicePlugin:
    """플러그인 팩토리 함수"""
    return SgNIPTPlugin()
