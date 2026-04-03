"""
Carrier screening: align fastq/analysis/output/log paths with CARRIER_SCREENING_FASTQ_DIR layout
and infer sequencing folder from FASTQ paths (vs portal sample_name / order id).
"""

from __future__ import annotations

import logging
import os
from typing import Optional, Tuple

from ...config import normalize_legacy_carrier_container_path, settings
from ...models import Job

logger = logging.getLogger(__name__)


def carrier_fastq_dir_from_local_paths(job: Job) -> Optional[str]:
    """
    When R1/R2 are local files in the same directory, use that directory as fastq_dir.

    This avoids ``work_root/fastq/<job.work_dir>/<sequencing_folder>`` when reusing FASTQ from
    an older run (different work segment under ``.../fastq/<work>/...``).
    """
    r1 = (job.fastq_r1_path or "").strip()
    r2 = (job.fastq_r2_path or "").strip()
    if not r1 or not r2:
        return None
    if not os.path.isfile(r1) or not os.path.isfile(r2):
        return None
    d1 = os.path.realpath(os.path.dirname(os.path.abspath(r1)))
    d2 = os.path.realpath(os.path.dirname(os.path.abspath(r2)))
    if d1 != d2:
        return None
    if not os.path.isdir(d1):
        return None
    return d1


def carrier_run_analysis_work_arg(job: Job) -> str:
    """
    ``-w`` for ``run_analysis.sh``: the wrapper expects ``fastq/<w>/<sequencing_folder>/``.

    When FASTQ lives under ``<work_root>/fastq/<work>/<seq>/``, use that ``<work>`` even if
    ``job.work_dir`` is a new date segment for a follow-up order.
    """
    wk = str(job.work_dir).strip() or "00"
    work_root = os.path.realpath(settings.carrier_screening_work_root)
    fastq_root = os.path.realpath(os.path.join(work_root, "fastq"))
    r1 = (job.fastq_r1_path or "").strip()
    if not r1 or not os.path.isfile(r1):
        return wk
    r1_dir = os.path.realpath(os.path.dirname(os.path.abspath(r1)))
    try:
        rel = os.path.relpath(r1_dir, fastq_root)
    except ValueError:
        return wk
    rel_norm = rel.replace("\\", "/")
    if rel_norm == "." or rel_norm.startswith("../"):
        return wk
    parts = [p for p in rel_norm.split("/") if p]
    if len(parts) >= 1:
        return parts[0]
    return wk


def carrier_sequencing_folder(job: Job) -> str:
    """
    fastq/<work_dir>/ 하위 디렉토리 이름 = run_analysis.sh --sample 인자.

    우선순위:
      1. params.carrier.sequencing_folder / fastq_sample_folder (명시적 지정)
      2. job.order_id (기본값)

    주의: R1 FASTQ 경로의 부모 디렉토리 이름으로 추론하지 않는다.
    FASTQ가 다른 order의 디렉토리에 있을 경우 잘못된 --sample이 전달되어
    다른 order의 FASTQ로 분석되고 결과가 덮어써지는 데이터 오염이 발생하기 때문.
    FASTQ는 반드시 fastq/{work_dir}/{order_id}/ 경로에 위치해야 한다.
    """
    carrier = (job.params or {}).get("carrier") or {}
    sf = carrier.get("sequencing_folder") or carrier.get("fastq_sample_folder")
    if isinstance(sf, str) and sf.strip():
        return sf.strip()
    return str(job.order_id)


def parse_carrier_output_tree_from_main_vcf(main_vcf: str) -> Tuple[str, str]:
    """
    main VCF가 .../output/<work_dir>/<leaf>/vcf/<file>.vcf.gz 에 있을 때
    (output_dir, tail) 를 반환. tail = '<work_dir>/<leaf>' (양 끝 슬래시 없음).

    그렇지 않으면 ("", "").
    """
    if not main_vcf:
        return "", ""
    main_vcf = normalize_legacy_carrier_container_path(main_vcf) or main_vcf
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
    main_vcf = normalize_legacy_carrier_container_path(main_vcf) or main_vcf
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
    # fastq / analysis / output / log 모두 order_id 기준으로 통일.
    # carrier_sequencing_folder 도 order_id 를 반환하므로 일관성 유지.
    seq_folder = carrier_sequencing_folder(job)   # == order_id (명시 재정의 없는 한)
    explicit_fq = carrier_fastq_dir_from_local_paths(job)
    new_f = explicit_fq if explicit_fq else os.path.join(work_root, "fastq", wk, seq_folder)
    new_a = os.path.join(work_root, "analysis", wk, seq_folder)
    new_o = os.path.join(work_root, "output", wk, seq_folder)
    new_l = os.path.join(work_root, "log", wk, seq_folder)
    changed = (
        job.fastq_dir != new_f
        or job.analysis_dir != new_a
        or job.output_dir != new_o
        or job.log_dir != new_l
    )
    if changed:
        logger.info(
            "[carrier_screening] Normalized paths for order %s (work=%s sample=%s) "
            "fastq=%s analysis=%s",
            job.order_id,
            wk,
            seq_folder,
            new_f,
            new_a,
        )
    job.fastq_dir = new_f
    job.analysis_dir = new_a
    job.output_dir = new_o
    job.log_dir = new_l
    return changed
