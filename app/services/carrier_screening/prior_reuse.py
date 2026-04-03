"""
Carrier: reuse completed pipeline outputs from a prior order (same sequencing run, new panel / order ID).

Activated when params.carrier.reuse_prior_pipeline_outputs is true and prior_order_id is set.
"""

from __future__ import annotations

import logging
import os
from typing import Optional, Tuple

from ...config import normalize_legacy_carrier_container_path
from ...models import Job, OrderStatus
from ...queue_manager import get_queue_manager

logger = logging.getLogger(__name__)


def carrier_reuse_prior_pipeline_requested(job: Job) -> bool:
    carrier = (job.params or {}).get("carrier") or {}
    return bool(carrier.get("reuse_prior_pipeline_outputs"))


def carrier_prior_order_id(job: Job) -> str:
    carrier = (job.params or {}).get("carrier") or {}
    pid = carrier.get("prior_order_id")
    if isinstance(pid, str) and pid.strip():
        return pid.strip()
    return ""


def resolve_carrier_prior_job(prior_order_id: str) -> Optional[Job]:
    """Memory first, then SQLite."""
    qm = get_queue_manager()
    j = qm.get_job(prior_order_id)
    if j:
        return j
    store = qm.store
    if store:
        return store.fetch_job(prior_order_id)
    return None


def _prior_has_pipeline_paths(prior: Job) -> bool:
    for key in ("analysis_dir", "output_dir"):
        p = getattr(prior, key, None) or (prior.model_dump().get(key) if prior else None)
        if isinstance(p, str) and p.strip() and os.path.isdir(p.strip()):
            return True
    return False


def validate_prior_for_pipeline_reuse(prior: Job) -> Tuple[bool, str]:
    # Same Nextflow stack as carrier_screening; whole_exome / health_screening store VCF/output the same way.
    if prior.service_code not in ("carrier_screening", "whole_exome", "health_screening"):
        return (
            False,
            f"prior order must be carrier_screening, whole_exome, or health_screening (got {prior.service_code!r})",
        )
    if prior.status not in (OrderStatus.COMPLETED, OrderStatus.REPORT_READY):
        return (
            False,
            f"prior order must be COMPLETED or REPORT_READY (got {prior.status.value})",
        )
    if not _prior_has_pipeline_paths(prior):
        return (
            False,
            "prior order has no analysis/output directories on disk — run the full pipeline first",
        )
    return True, ""


def apply_carrier_prior_reuse_metadata(job: Job, prior: Job) -> None:
    """Stash paths for VCF search, QC, and pipeline_complete fallback."""
    job.params["_prior_reuse"] = True
    job.params["_prior_reuse_analysis_dir"] = prior.analysis_dir
    job.params["_prior_reuse_output_dir"] = prior.output_dir
    job.params["_prior_reuse_log_dir"] = prior.log_dir
    job.params["_prior_reuse_sample_name"] = prior.sample_name
    job.params["_prior_reuse_order_id"] = prior.order_id
    main = (prior.params or {}).get("main_vcf")
    if isinstance(main, str) and main.strip():
        main = normalize_legacy_carrier_container_path(main.strip()) or main.strip()
        if os.path.isfile(main):
            job.params["_prior_reuse_main_vcf_hint"] = os.path.abspath(main)


def prepare_carrier_prior_pipeline_reuse(job: Job) -> bool:
    """
    If reuse is requested, validate prior order and set job.params metadata.
    Returns True when the caller should skip FASTQ prep; False when reuse is not requested.
    Raises RuntimeError when reuse is requested but cannot be satisfied (caller should surface message).
    """
    if not carrier_reuse_prior_pipeline_requested(job):
        return False
    pid = carrier_prior_order_id(job)
    if not pid:
        msg = "reuse_prior_pipeline_outputs is set but prior_order_id is empty"
        logger.error("[carrier_screening] %s", msg)
        raise RuntimeError(msg)
    prior = resolve_carrier_prior_job(pid)
    if not prior:
        msg = f"prior order not found: {pid!r} (save or run the prior order so it exists in the daemon DB)"
        logger.error("[carrier_screening] %s", msg)
        raise RuntimeError(msg)
    ok, err = validate_prior_for_pipeline_reuse(prior)
    if not ok:
        logger.error(
            "[carrier_screening] prior order %s invalid for reuse: %s", pid, err
        )
        raise RuntimeError(f"Prior pipeline reuse: {err}")

    apply_carrier_prior_reuse_metadata(job, prior)
    logger.info(
        "[carrier_screening] prior pipeline reuse: new=%s prior=%s (analysis=%s)",
        job.order_id,
        pid,
        prior.analysis_dir,
    )
    return True
