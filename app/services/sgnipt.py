"""
sgNIPT Service Plugin

Docker 이미지 `sgnipt`에서 run_sgnipt.sh 를 실행합니다.
호스트 레이아웃: {job_root}/fastq|analysis|log|output/<work_dir>/<order_id>/
(job_root = SGNIPT_WORK_ROOT 또는 SGNIPT_LAYOUT_ROOT)
완료 판정: output/.../<order_id>/<order_id>.json (generate_result_json.py)
Portal 연동: 동일 JSON을 output/.../result.json 으로 복사합니다.
"""

import os
import glob
import json
import logging
import shlex
import shutil
from typing import Any, Dict, List, Tuple

from .base import ServicePlugin
from app.config import settings
from app.models import Job, OutputFile

logger = logging.getLogger(__name__)


def apply_sgnipt_layout_directories(job: Job) -> bool:
    """
    Job의 fastq_dir / analysis_dir / output_dir / log_dir 를 현재 설정 기준으로 갱신합니다.
    DB에 저장된 구 경로를 현재 SGNIPT_WORK_ROOT 로 교체할 때 사용.
    변경이 있었으면 True 반환.
    """
    if job.service_code != "sgnipt":
        return False
    root = settings.sgnipt_job_root
    wk = str(job.work_dir or "").strip() or "00"
    oid = str(job.order_id or "").strip()
    new_f = os.path.join(root, "fastq", wk, oid)
    new_a = os.path.join(root, "analysis", wk, oid)
    new_o = os.path.join(root, "output", wk, oid)
    new_l = os.path.join(root, "log", wk, oid)
    changed = (
        job.fastq_dir != new_f
        or job.analysis_dir != new_a
        or job.output_dir != new_o
        or job.log_dir != new_l
    )
    if changed:
        logger.info(
            "[sgnipt] Normalized paths for %s (work=%s root=%s)",
            job.order_id, wk, root,
        )
        job.fastq_dir = new_f
        job.analysis_dir = new_a
        job.output_dir = new_o
        job.log_dir = new_l
    return changed


class SgNIPTPlugin(ServicePlugin):
    """sgNIPT 파이프라인 (Docker 또는 호스트 bash)."""

    @property
    def service_code(self) -> str:
        return "sgnipt"

    @property
    def display_name(self) -> str:
        return "sgNIPT (Single-Gene NIPT)"

    def validate_params(self, params: Dict[str, Any]) -> Tuple[bool, str]:
        return True, ""

    def _order_result_json(self, job: Job) -> str:
        oid = (job.order_id or "").strip()
        return os.path.join(job.output_dir, f"{oid}.json")

    async def prepare_inputs(self, job: Job) -> bool:
        try:
            for d in (job.fastq_dir, job.analysis_dir, job.output_dir, job.log_dir):
                if d:
                    os.makedirs(d, exist_ok=True)

            if job.fastq_r1_path and job.fastq_r2_path:
                for src in [job.fastq_r1_path, job.fastq_r2_path]:
                    dst = os.path.join(job.fastq_dir, os.path.basename(src))
                    if not os.path.exists(dst):
                        os.symlink(os.path.abspath(src), dst)

            logger.info("[sgnipt] Input prepared for %s (work_dir=%s)", job.order_id, job.work_dir)
            return True
        except OSError as e:
            logger.error("[sgnipt] Failed to prepare inputs: %s", e)
            return False

    def _pipeline_command_parts(self, job: Job) -> List[str]:
        script = settings.sgnipt_run_script_path
        oid = (job.order_id or "").strip()
        wd = (job.work_dir or "").strip()
        # run_sgnipt.sh itself calls docker run (like carrier's run_analysis.sh).
        # Call it directly from the daemon — docker socket is already mounted.
        return ["bash", script, "--order_id", oid, "--work_dir", wd]

    async def get_pipeline_command(self, job: Job) -> str:
        # SGNIPT_ROOT_DIR must be available as an env var when the command runs
        # (set in .env.docker and inherited by the subprocess).
        parts = self._pipeline_command_parts(job)
        return shlex.join(parts)

    async def check_completion(self, job: Job) -> bool:
        path = self._order_result_json(job)
        if not os.path.isfile(path):
            logger.warning("[sgnipt] Completion file missing: %s", path)
            return False
        try:
            with open(path, "r", encoding="utf-8") as f:
                json.load(f)
        except (OSError, json.JSONDecodeError) as e:
            logger.warning("[sgnipt] Invalid or unreadable result JSON %s: %s", path, e)
            return False
        return True

    async def process_results(self, job: Job) -> bool:
        src = self._order_result_json(job)
        dst = os.path.join(job.output_dir, "result.json")
        try:
            if not os.path.isfile(src):
                logger.error("[sgnipt] Expected result not found: %s", src)
                return False
            if os.path.exists(dst) and os.path.samefile(src, dst):
                logger.info("[sgnipt] result.json already same as %s", os.path.basename(src))
                return True
            shutil.copy2(src, dst)
            logger.info("[sgnipt] Copied order result to result.json for portal")
            return True
        except OSError as e:
            logger.error("[sgnipt] process_results failed: %s", e)
            return False

    async def get_output_files(self, job: Job) -> List[OutputFile]:
        out = job.output_dir
        files: List[OutputFile] = []
        oid = (job.order_id or "").strip()

        result_json = os.path.join(out, "result.json")
        if os.path.isfile(result_json):
            files.append(
                OutputFile(
                    file_path=result_json,
                    file_type="review_json",
                    file_name="result.json",
                    content_type="application/json",
                )
            )

        order_json = os.path.join(out, f"{oid}.json")
        if os.path.isfile(order_json) and os.path.abspath(order_json) != os.path.abspath(result_json):
            files.append(
                OutputFile(
                    file_path=order_json,
                    file_type="order_result_json",
                    file_name=f"{oid}.json",
                    content_type="application/json",
                )
            )

        tar_path = os.path.join(out, f"{oid}.output.tar")
        if os.path.isfile(tar_path):
            files.append(
                OutputFile(
                    file_path=tar_path,
                    file_type="output_tar",
                    file_name=f"{oid}.output.tar",
                    content_type="application/x-tar",
                )
            )

        mq = os.path.join(out, "multiqc", "multiqc_report.html")
        if os.path.isfile(mq):
            files.append(
                OutputFile(
                    file_path=mq,
                    file_type="multiqc_html",
                    file_name="multiqc_report.html",
                    content_type="text/html",
                )
            )

        for f in glob.glob(os.path.join(out, "**", "*_final_report.html"), recursive=True):
            if os.path.isfile(f):
                files.append(
                    OutputFile(
                        file_path=f,
                        file_type="final_report_html",
                        file_name=os.path.basename(f),
                        content_type="text/html",
                    )
                )

        logger.info("[sgnipt] Output files for upload: %s", len(files))
        return files

    def get_progress_stages(self) -> Dict[str, int]:
        return {
            "FASTQ": 15,
            "ALIGN": 35,
            "VARIANT": 60,
            "REPORT": 85,
            "SUMMARY": 92,
        }


def create_plugin() -> ServicePlugin:
    return SgNIPTPlugin()
