"""
Carrier screening: align fastq/analysis/output/log paths with CARRIER_SCREENING_FASTQ_DIR layout
and infer sequencing folder from FASTQ paths (vs portal sample_name / order id).
"""

from __future__ import annotations

import logging
import os
from ...config import settings
from ...models import Job

logger = logging.getLogger(__name__)


def carrier_sequencing_folder(job: Job) -> str:
    """
    Directory name under fastq/<work_dir>/ that matches on-disk FASTQ layout.
    Prefer explicit params; else parent dir of R1; else job.sample_name.
    """
    carrier = (job.params or {}).get("carrier") or {}
    sf = carrier.get("sequencing_folder") or carrier.get("fastq_sample_folder")
    if isinstance(sf, str) and sf.strip():
        return sf.strip()
    r1 = job.fastq_r1_path or ""
    if r1 and os.path.isfile(r1):
        folder = os.path.basename(os.path.dirname(os.path.abspath(r1)))
        if folder and folder != str(job.sample_name):
            logger.info(
                "[carrier_screening] sequencing folder=%s (from FASTQ path; job.sample_name=%s)",
                folder,
                job.sample_name,
            )
        if folder:
            return folder
    return str(job.sample_name)


def apply_carrier_layout_directories(job: Job) -> bool:
    """
    Set fastq paths under carrier_screening_layout_base; analysis/output/log under
    carrier_screening_work_root (CARRIER_SCREENING_ARTIFACT_BASE or layout base).

    - fastq_dir: layout_base/fastq/<work>/<sequencing_folder>/
    - analysis|output|log: work_root/<kind>/<work>/<sample_name>/

    Returns True if any path field changed.
    """
    if job.service_code != "carrier_screening":
        return False
    layout_base = settings.carrier_screening_layout_base
    work_root = settings.carrier_screening_work_root
    wk = str(job.work_dir).strip() or "00"
    seq_folder = carrier_sequencing_folder(job)
    path_key = str(job.sample_name).strip() or job.order_id
    new_f = os.path.join(layout_base, "fastq", wk, seq_folder)
    new_a = os.path.join(work_root, "analysis", wk, path_key)
    new_o = os.path.join(work_root, "output", wk, path_key)
    new_l = os.path.join(work_root, "log", wk, path_key)
    changed = (
        job.fastq_dir != new_f
        or job.analysis_dir != new_a
        or job.output_dir != new_o
        or job.log_dir != new_l
    )
    if changed:
        logger.info(
            "[carrier_screening] Normalized paths for order %s (work=%s seq=%s sample=%s) → analysis=%s",
            job.order_id,
            wk,
            seq_folder,
            path_key,
            new_a,
        )
    job.fastq_dir = new_f
    job.analysis_dir = new_a
    job.output_dir = new_o
    job.log_dir = new_l
    return changed
