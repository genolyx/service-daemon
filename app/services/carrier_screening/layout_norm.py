"""
Carrier screening: align fastq/analysis/output/log paths with CARRIER_SCREENING_FASTQ_DIR layout
and infer sequencing folder from FASTQ paths (vs portal sample_name / order id).
"""

from __future__ import annotations

import logging
import os
from typing import Tuple

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


def parse_carrier_output_tree_from_main_vcf(main_vcf: str) -> Tuple[str, str]:
    """
    main VCF가 .../output/<work_dir>/<leaf>/vcf/<file>.vcf.gz 에 있을 때
    (output_dir, tail) 를 반환. tail = '<work_dir>/<leaf>' (양 끝 슬래시 없음).

    그렇지 않으면 ("", "").
    """
    if not main_vcf:
        return "", ""
    m = os.path.realpath(os.path.abspath(main_vcf))
    vdir = os.path.dirname(m)
    if os.path.basename(vdir) != "vcf":
        return "", ""
    output_dir = os.path.dirname(vdir)
    rel = output_dir.replace("\\", "/")
    marker = "/output/"
    if marker not in rel:
        return "", ""
    prefix, _, tail = rel.partition(marker)
    tail = (tail or "").strip("/")
    if not prefix or not tail:
        return "", ""
    return output_dir, tail


def align_carrier_job_dirs_from_main_vcf(job: Job, main_vcf: str) -> bool:
    """
    check_completion 이 layout_base·script_data 쪽에서 main_vcf 를 찾았는데
    job.output_dir 은 work_root 기준으로 다른 트리를 가리키는 경우(흔한 Docker/이중 루트)에
    main VCF가 있는 실제 output/ 분석/ log 트리로 analysis|output|log 를 맞춘다.

    Returns:
        경로가 하나라도 바뀌었으면 True.
    """
    if job.service_code != "carrier_screening":
        return False
    output_dir, tail = parse_carrier_output_tree_from_main_vcf(main_vcf)
    if not output_dir or not tail:
        return False
    rel = output_dir.replace("\\", "/")
    marker = "/output/"
    prefix, _, _ = rel.partition(marker)
    analysis_dir = f"{prefix}/analysis/{tail}"
    log_dir = f"{prefix}/log/{tail}"

    changed = False
    if job.output_dir != output_dir:
        job.output_dir = output_dir
        changed = True
    if os.path.isdir(analysis_dir) and job.analysis_dir != analysis_dir:
        job.analysis_dir = analysis_dir
        changed = True
    if os.path.isdir(log_dir) and job.log_dir != log_dir:
        job.log_dir = log_dir
        changed = True
    if changed:
        logger.info(
            "[carrier_screening] Artifact dirs aligned from main VCF for order %s → output=%s analysis=%s",
            job.order_id,
            output_dir,
            job.analysis_dir,
        )
    return changed


def apply_carrier_layout_directories(job: Job) -> bool:
    """
    Set fastq and analysis/output/log under carrier_screening_work_root
    (CARRIER_SCREENING_ARTIFACT_BASE when set, else layout base).

    - fastq_dir: work_root/fastq/<work>/<sequencing_folder>/
    - analysis|output|log: work_root/<kind>/<work>/<sample_name>/

    Returns True if any path field changed.
    """
    if job.service_code != "carrier_screening":
        return False
    work_root = settings.carrier_screening_work_root
    wk = str(job.work_dir).strip() or "00"
    seq_folder = carrier_sequencing_folder(job)
    path_key = str(job.sample_name).strip() or job.order_id
    new_f = os.path.join(work_root, "fastq", wk, seq_folder)
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
