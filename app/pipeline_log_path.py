"""
Resolve on-disk pipeline log files for Portal ``GET /order/{id}/pipeline-log``.

Carrier / whole exome / health screening (gx-exome):
  ``{layout_base}/log/<work>/<sample>/nextflow.log``

sgNIPT:
  ``{job_root}/log/<work>/<order_id>/nextflow.log`` (layout ↔ work mount aliases)

Other services (e.g. NIPT):
  ``{log_dir}/pipeline.log`` (legacy)
"""

from __future__ import annotations

import os
from typing import List, Optional, Tuple

from .config import settings
from .models import Job
from .services.carrier_screening.layout_norm import (
    _CARRIER_LIKE,
    carrier_run_analysis_work_arg,
    carrier_sequencing_folder,
)

_GX_EXOME_LAYOUT_PREFIX = "/data/gx-exome"
_LOG_NAME_NEXTFLOW = "nextflow.log"
_LOG_NAME_PIPELINE = "pipeline.log"


def _add_candidate(seen: set[str], out: List[str], path: str) -> None:
    p = (path or "").strip()
    if not p:
        return
    try:
        ap = os.path.abspath(p)
    except OSError:
        return
    if ap in seen:
        return
    seen.add(ap)
    out.append(ap)


def _carrier_log_candidates(job: Job) -> List[str]:
    work = carrier_run_analysis_work_arg(job)
    sample = carrier_sequencing_folder(job)
    rel = os.path.join("log", work, sample, _LOG_NAME_NEXTFLOW)
    seen: set[str] = set()
    out: List[str] = []

    layout = (settings.carrier_screening_layout_base or "").strip().rstrip("/")
    if layout:
        _add_candidate(seen, out, os.path.join(layout, rel))

    script_data = (settings.carrier_screening_script_data_dir or "").strip().rstrip("/")
    if script_data:
        _add_candidate(seen, out, os.path.join(script_data, rel))

    host = (settings.carrier_screening_host or "").strip().rstrip("/")
    if host:
        _add_candidate(seen, out, os.path.join(host, rel))
        if layout.startswith(_GX_EXOME_LAYOUT_PREFIX):
            suffix = rel
            if host and layout == _GX_EXOME_LAYOUT_PREFIX:
                alias = os.path.join(_GX_EXOME_LAYOUT_PREFIX, suffix)
                _add_candidate(seen, out, alias)

    log_dir = (job.log_dir or "").strip()
    if log_dir:
        _add_candidate(seen, out, os.path.join(log_dir, _LOG_NAME_NEXTFLOW))

    return out


def _sgnipt_log_candidates(job: Job) -> List[str]:
    wk = str(job.work_dir or "").strip() or "00"
    oid = str(job.order_id or "").strip()
    seen: set[str] = set()
    out: List[str] = []

    for base in (
        (job.log_dir or "").strip(),
        os.path.join(settings.sgnipt_job_root, "log", wk, oid),
        os.path.join((settings.sgnipt_layout_root or "").strip(), "log", wk, oid),
    ):
        if base:
            _add_candidate(seen, out, os.path.join(base, _LOG_NAME_NEXTFLOW))

    return out


def order_pipeline_log_candidates(job: Job) -> List[Tuple[str, str]]:
    """
    Ordered (absolute_path, display_name) pairs to probe.
    Only ``nextflow.log`` for gx-exome / sgNIPT (not ``nextflow.log.1`` rotations).
    """
    svc = (job.service_code or "").strip()
    if svc in _CARRIER_LIKE:
        return [(p, _LOG_NAME_NEXTFLOW) for p in _carrier_log_candidates(job)]
    if svc == "sgnipt":
        return [(p, _LOG_NAME_NEXTFLOW) for p in _sgnipt_log_candidates(job)]
    log_dir = (job.log_dir or "").strip()
    if log_dir:
        return [(os.path.join(log_dir, _LOG_NAME_PIPELINE), _LOG_NAME_PIPELINE)]
    return []


def resolve_order_pipeline_log(job: Job) -> Tuple[Optional[str], str]:
    """Return (existing file path, log filename for UI/errors)."""
    candidates = order_pipeline_log_candidates(job)
    if not candidates:
        return None, _LOG_NAME_PIPELINE
    display = candidates[0][1]
    for path, name in candidates:
        if os.path.isfile(path):
            return path, name
    return None, display
