"""주문 런 삭제: analysis / output / log 트리만 제거 (FASTQ 디렉터리는 건드리지 않음)."""

from __future__ import annotations

import logging
import os
import shutil
from typing import List, Set, Tuple

from .config import settings
from .models import Job
from .services.carrier_screening.layout_norm import carrier_sequencing_folder

logger = logging.getLogger(__name__)


def collect_run_artifact_directories(job: Job) -> List[str]:
    """Job에 연결된 분석 산출물 경로 후보 (존재하는 디렉터리만)."""
    dirs: List[str] = []
    for attr in ("analysis_dir", "output_dir", "log_dir"):
        p = getattr(job, attr, None)
        if isinstance(p, str) and p.strip():
            dirs.append(os.path.abspath(p.strip()))

    if job.service_code == "carrier_screening":
        lb = settings.carrier_screening_layout_base
        wk = str(job.work_dir).strip() or "00"
        sk = str(job.sample_name).strip() or job.order_id
        seq = carrier_sequencing_folder(job)
        folder_keys = list(dict.fromkeys([x for x in (sk, seq) if x]))
        extra_roots: List[str] = [lb]
        sd = (settings.carrier_screening_script_data_dir or "").strip()
        if sd:
            extra_roots.append(os.path.abspath(sd))
        for root in extra_roots:
            for fk in folder_keys:
                for sub in ("analysis", "output", "log"):
                    ap = os.path.abspath(os.path.join(root, sub, wk, fk))
                    if ap not in dirs:
                        dirs.append(ap)

    seen_real: Set[str] = set()
    out: List[str] = []
    for d in dirs:
        if not os.path.isdir(d):
            continue
        try:
            r = os.path.realpath(d)
        except OSError:
            continue
        if r in seen_real:
            continue
        seen_real.add(r)
        out.append(d)
    return out


def _allowed_roots() -> List[str]:
    roots: List[str] = []
    for val in (
        settings.carrier_screening_work_root,
        settings.carrier_screening_layout_base,
        settings.base_dir,
        settings.analysis_base_dir,
        settings.output_base_dir,
        settings.log_base_dir,
        settings.carrier_screening_script_data_dir,
    ):
        if isinstance(val, str) and val.strip():
            try:
                roots.append(os.path.realpath(val.strip()))
            except OSError:
                continue
    return list(dict.fromkeys(roots))


def delete_run_artifacts(job: Job) -> Tuple[List[str], List[str]]:
    """
    Returns:
        (deleted_paths, error_messages)
    """
    deleted: List[str] = []
    errors: List[str] = []
    fq_real: str | None = None
    if job.fastq_dir and str(job.fastq_dir).strip():
        try:
            fq_real = os.path.realpath(os.path.abspath(str(job.fastq_dir).strip()))
        except OSError:
            fq_real = None

    allowed = _allowed_roots()
    segment_ok = ("analysis", "output", "log")

    for path in collect_run_artifact_directories(job):
        try:
            rp = os.path.realpath(path)
        except OSError as e:
            errors.append(f"{path}: {e}")
            continue
        if fq_real and (rp == fq_real or rp.startswith(fq_real + os.sep)):
            logger.warning("Skip delete (under fastq_dir): %s", path)
            continue
        if not any(rp == ar or rp.startswith(ar + os.sep) for ar in allowed):
            errors.append(f"{path}: refused (outside configured daemon directories)")
            continue
        parts = rp.split(os.sep)
        if not any(seg in parts for seg in segment_ok):
            errors.append(f"{path}: refused (path must include analysis, output, or log)")
            continue
        try:
            shutil.rmtree(path)
            deleted.append(path)
            logger.info("Deleted run artifact tree: %s", path)
        except OSError as e:
            errors.append(f"{path}: {e}")
    return deleted, errors
