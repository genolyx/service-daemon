"""
sgNIPT Service Plugin

Docker 이미지 `sgnipt`에서 run_sgnipt.sh 를 실행합니다.
호스트 레이아웃: {job_root}/fastq|analysis|log|output/<work_dir>/<order_id>/
(job_root = SGNIPT_WORK_ROOT; FASTQ 는 레이아웃 sgNIPT/fastq 와 동일 트리 — compose 가
sgnipt_work/fastq 에 레이아웃 fastq 를 마운트. run_sgnipt.sh 는 기본 포그라운드 docker run
— carrier run_analysis.sh 와 동일하게 쉘이 파이프라인 종료까지 블록함)
완료 판정: output/.../<order_id>/<order_id>.json (generate_result_json.py)
Portal 연동: 동일 JSON을 output/.../result.json 으로 복사합니다.
"""

import os
import glob
import json
import logging
import shlex
import shutil
from typing import Any, Dict, List, Optional, Set, Tuple

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
        return "Single-gene NIPT"

    def validate_params(self, params: Dict[str, Any], strict: bool = True) -> Tuple[bool, str]:
        return True, ""

    def sync_is_complete(self, job: Job) -> bool:
        """
        daemon 재시작 복구 시 파이프라인이 실제로 완료됐는지 동기적으로 확인.
        output_dir/<order_id>.json 파이프라인 결과 파일이 존재하면 완료로 판정.
        """
        path = self._order_result_json(job)
        if os.path.isfile(path):
            logger.info(
                "[sgnipt] sync_is_complete: found %s → marking COMPLETED", path
            )
            return True
        return False

    def _run_script_candidates(self) -> List[str]:
        """
        존재 여부와 무관한 후보 목록 (오류 메시지·탐색용).
        Docker 에서 소스(/home/ken/sgnipt) 와 데이터 트리(/home/ken/sgNIPT) 가 갈라질 수 있음.
        """
        candidates: List[str] = []
        configured = (settings.sgnipt_run_script_path or "").strip()
        if configured:
            candidates.append(os.path.abspath(configured))
        src_root = (settings.sgnipt_src_root or "").strip()
        if src_root:
            candidates.append(
                os.path.abspath(os.path.join(src_root, "src", "run_sgnipt.sh"))
            )
        root = (settings.sgnipt_job_root or "").strip()
        if root:
            candidates.append(os.path.abspath(os.path.join(root, "src", "run_sgnipt.sh")))
        layout = (settings.sgnipt_layout_root or "").strip()
        if layout:
            candidates.append(os.path.abspath(os.path.join(layout, "src", "run_sgnipt.sh")))
        seen: Set[str] = set()
        ordered: List[str] = []
        for c in candidates:
            if not c or c in seen:
                continue
            seen.add(c)
            ordered.append(c)
        return ordered

    def _resolve_run_script(self) -> str:
        """
        SGNIPT_RUN_SCRIPT_PATH 가 틀리거나(대소문자 경로 등) 없을 때 job / layout / SGNIPT_SRC_ROOT 로 보완.
        """
        configured = (settings.sgnipt_run_script_path or "").strip()
        for c in self._run_script_candidates():
            if os.path.isfile(c):
                if configured and os.path.abspath(configured) != c:
                    logger.warning(
                        "[sgnipt] Using run script %s (configured path missing: %s)",
                        c,
                        configured,
                    )
                return c
        cand = self._run_script_candidates()
        return cand[0] if cand else ""

    def get_pipeline_cwd(self, job: Job) -> Optional[str]:
        """
        수동 실행과 동일하게 **소스 저장소 루트**에서 ``bash src/run_sgnipt.sh ...`` 하도록 cwd 고정.
        SGNIPT_WORK_ROOT 만 쓰면 data 트리에 src/ 가 없어 상대 경로가 깨질 수 있음.
        """
        script = self._resolve_run_script()
        if script and os.path.isfile(script):
            repo_root = os.path.abspath(os.path.join(os.path.dirname(script), ".."))
            if os.path.isdir(os.path.join(repo_root, "src")):
                return repo_root
        root = (settings.sgnipt_job_root or "").strip()
        if root and os.path.isdir(root):
            return os.path.abspath(root)
        return None

    def _order_result_json(self, job: Job) -> str:
        oid = (job.order_id or "").strip()
        return os.path.join(job.output_dir, f"{oid}.json")

    def _fastq_host_root(self) -> str:
        """run_sgnipt.sh 의 HOST_FASTQ_DIR (= SGNIPT_ROOT_DIR/fastq) 과 동일해야 함."""
        root = (settings.sgnipt_job_root or "").strip()
        return os.path.abspath(os.path.join(root, "fastq"))

    @staticmethod
    def _file_is_under_dir(dir_abs: str, path_abs: str) -> bool:
        dir_abs = os.path.abspath(dir_abs)
        path_abs = os.path.abspath(path_abs)
        try:
            return os.path.commonpath([dir_abs, path_abs]) == dir_abs
        except ValueError:
            return False

    def _run_sgnipt_fastq_rel_paths(self, job: Job) -> Tuple[Optional[str], Optional[str]]:
        """
        run_sgnipt / inner Docker 가 기대하는 FASTQ 경로: HOST_FASTQ_DIR 기준 상대 경로.
        """
        r1 = (job.fastq_r1_path or "").strip()
        r2 = (job.fastq_r2_path or "").strip()
        if not r1 or not r2:
            return None, None
        fastq_host = self._fastq_host_root()
        dest1 = os.path.abspath(os.path.join(job.fastq_dir or "", os.path.basename(r1)))
        dest2 = os.path.abspath(os.path.join(job.fastq_dir or "", os.path.basename(r2)))
        if not os.path.isfile(dest1) or not os.path.isfile(dest2):
            return None, None
        try:
            rel1 = os.path.relpath(dest1, fastq_host)
            rel2 = os.path.relpath(dest2, fastq_host)
        except ValueError:
            return None, None
        if rel1.startswith("..") or rel2.startswith(".."):
            return None, None
        return rel1.replace("\\", "/"), rel2.replace("\\", "/")

    async def prepare_inputs(self, job: Job) -> bool:
        for d in (job.fastq_dir, job.analysis_dir, job.output_dir, job.log_dir):
            if d:
                try:
                    os.makedirs(d, exist_ok=True)
                except OSError as e:
                    raise RuntimeError(
                        f"sgnipt: cannot create job directory {d}: {e}"
                    ) from e

        input_bam_csv = (job.params or {}).get("input_bam_csv", "").strip()
        if input_bam_csv:
            if not os.path.isfile(input_bam_csv):
                raise RuntimeError(
                    f"sgnipt: BAM samplesheet CSV not found at path visible to daemon: {input_bam_csv}"
                )
            logger.info("[sgnipt] BAM simulation mode — using CSV: %s", input_bam_csv)
        else:
            r1 = (job.fastq_r1_path or "").strip()
            r2 = (job.fastq_r2_path or "").strip()
            fqdir = os.path.abspath(job.fastq_dir or "")
            if r1 and r2:
                for label, src in (("R1", r1), ("R2", r2)):
                    src_abs = os.path.abspath(src)
                    if not os.path.isfile(src_abs):
                        raise RuntimeError(
                            f"sgnipt: {label} FASTQ not found at path visible to daemon: {src_abs}"
                        )
                    if fqdir and self._file_is_under_dir(fqdir, src_abs):
                        continue
                    dst = os.path.join(job.fastq_dir or "", os.path.basename(src_abs))
                    if os.path.lexists(dst):
                        continue
                    try:
                        os.symlink(src_abs, dst)
                    except OSError as e:
                        raise RuntimeError(
                            f"sgnipt: could not symlink FASTQ into {dst}: {e}"
                        ) from e
            elif r1 or r2:
                raise RuntimeError(
                    "sgnipt: both fastq_r1_path and fastq_r2_path are required"
                )

        script = self._resolve_run_script()
        if not script or not os.path.isfile(script):
            tried = ", ".join(repr(p) for p in self._run_script_candidates())
            raise RuntimeError(
                "sgnipt: run_sgnipt.sh not found. For Docker, set SGNIPT_SRC_ROOT to the "
                f"source clone (same path as SGNIPT_SRC_HOST), or set SGNIPT_RUN_SCRIPT_PATH. Tried: {tried}"
            )

        logger.info("[sgnipt] Input prepared for %s (work_dir=%s)", job.order_id, job.work_dir)
        return True

    def _pipeline_command_parts(self, job: Job) -> List[str]:
        raw = self._resolve_run_script()
        script = os.path.abspath(raw) if raw else ""
        oid = (job.order_id or "").strip()
        wd = (job.work_dir or "").strip()
        params = job.params or {}

        input_bam_csv = params.get("input_bam_csv", "").strip()
        if input_bam_csv:
            # BAM simulation mode: pass CSV path directly to run_sgnipt.sh --input-bam
            parts: List[str] = [
                "bash", script,
                "--order-id", oid,
                "--work-id", wd,
                "--input-bam", input_bam_csv,
            ]
            logger.info("[sgnipt] BAM simulation command — input-bam=%s", input_bam_csv)
        else:
            # run_sgnipt.sh itself calls docker run (like carrier's run_analysis.sh).
            # cwd 는 get_pipeline_cwd → 저장소 루트(수동 실행과 동일).
            parts = ["bash", script, "--order_id", oid, "--work_dir", wd]
            fr1, fr2 = self._run_sgnipt_fastq_rel_paths(job)
            if fr1 and fr2:
                parts.extend(["--fastq_r1", fr1, "--fastq_r2", fr2])

        if params.get("_pipeline_fresh"):
            parts.append("--fresh")
        return parts

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

    @staticmethod
    def _find_variant_report_json(analysis_dir: Optional[str], output_dir: Optional[str], sample_id: str) -> Optional[str]:
        """
        variant_report.json 위치를 탐색합니다.
        파이프라인 버전/레이아웃에 따라 경로가 다를 수 있으므로
        고정 경로를 먼저 시도하고, 없으면 glob으로 탐색 (work/ 제외).
        """
        fixed = [
            os.path.join(analysis_dir or "", sample_id, "variants", f"{sample_id}.variant_report.json"),
            os.path.join(analysis_dir or "", "variant", f"{sample_id}.variant_report.json"),
            os.path.join(analysis_dir or "", "variants", f"{sample_id}.variant_report.json"),
            os.path.join(output_dir or "", sample_id, "variants", f"{sample_id}.variant_report.json"),
        ]
        found = next((p for p in fixed if p and os.path.isfile(p)), None)
        if found:
            return found
        for root in filter(None, [analysis_dir, output_dir]):
            hits = [
                p for p in glob.glob(os.path.join(root, "**", f"{sample_id}.variant_report.json"), recursive=True)
                if os.sep + "work" + os.sep not in p
            ]
            if hits:
                return sorted(hits)[0]
        return None

    @staticmethod
    def _build_merged_result(summary: dict, vr: dict) -> dict:
        """
        파이프라인 요약 JSON(summary)과 variant_report JSON(vr)을 병합하여
        포털 Review에서 바로 사용할 수 있는 완성된 result.json 구조를 반환합니다.
        carrier_screening의 result.json과 동일한 단일-파일 구조가 목표입니다.
        """
        samples = summary.get("samples") or []
        sample0 = samples[0] if samples else {}
        sample_id = vr.get("sample_id") or sample0.get("sample_id") or ""
        return {
            **summary,
            # variant_report 상세 필드 (Review 핵심 데이터)
            "sample_id": sample_id,
            "panel": vr.get("panel"),
            "fetal_fraction_used": vr.get("fetal_fraction_used"),
            "summary": vr.get("summary"),
            "clinical_findings": vr.get("clinical_findings") or [],
            "all_target_variants": vr.get("all_target_variants") or [],
            "gene_coverage_validation": vr.get("gene_coverage_validation"),
            # QC: samples[0] 하위 필드를 최상위로 올림
            "sgnipt_status": summary.get("status"),
            "sgnipt_status_flags": sample0.get("status_flags") or [],
            "fastq_qc": sample0.get("fastq_qc"),
            "bam_qc": sample0.get("bam_qc"),
            "fetal_fraction_detail": sample0.get("fetal_fraction"),
            "variant_analysis_summary": sample0.get("variant_analysis"),
        }

    async def process_results(self, job: Job) -> bool:
        src = self._order_result_json(job)
        dst = os.path.join(job.output_dir, "result.json")
        try:
            if not os.path.isfile(src):
                logger.error("[sgnipt] Expected result not found: %s", src)
                return False

            with open(src, "r", encoding="utf-8") as f:
                summary: dict = json.load(f)

            samples = summary.get("samples") or []
            sample0 = samples[0] if samples else {}
            sample_id = sample0.get("sample_id") or (job.order_id or "").strip()

            vr_path = self._find_variant_report_json(job.analysis_dir, job.output_dir, sample_id)
            vr: dict = {}
            if vr_path:
                logger.info("[sgnipt] variant_report.json found at %s", vr_path)
                try:
                    with open(vr_path, "r", encoding="utf-8") as f:
                        vr = json.load(f) or {}
                except Exception as e:
                    logger.warning("[sgnipt] Could not read variant_report.json %s: %s", vr_path, e)
            else:
                logger.warning(
                    "[sgnipt] variant_report.json not found (sample_id=%s, analysis_dir=%s)",
                    sample_id, job.analysis_dir,
                )

            merged = self._build_merged_result(summary, vr)
            text = json.dumps(merged, ensure_ascii=False, indent=2, default=str)
            with open(dst, "w", encoding="utf-8") as f:
                f.write(text)
            logger.info("[sgnipt] Merged result.json written to %s", dst)
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
