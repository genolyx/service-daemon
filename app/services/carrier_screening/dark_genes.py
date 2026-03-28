"""
Collect dark-gene pipeline outputs (unified Nextflow outdir = analysis_dir).

Looks for summary reports under summary/ or results/summary/ (dark_gene_pipeline layout).
"""

from __future__ import annotations

import glob
import logging
import os
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

_MAX_READ_BYTES = 500_000
_SUMMARY_GLOBS = (
    "summary/*_summary_report.txt",
    "summary/*summary_report*.txt",
    "results/summary/*_summary_report.txt",
    "results/summary/*summary_report*.txt",
)
_DETAILED_GLOBS = (
    "summary/*_detailed_report.txt",
    "summary/*detailed_report*.txt",
    "results/summary/*_detailed_report.txt",
)


def _read_text_safe(path: str, max_chars: int) -> str:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return f.read(max_chars + 1)[:max_chars]
    except OSError as e:
        logger.warning("[dark_genes] Could not read %s: %s", path, e)
        return ""


def _first_existing_glob(base: str, patterns: tuple) -> List[str]:
    out: List[str] = []
    seen = set()
    for pat in patterns:
        for p in glob.glob(os.path.join(base, pat)):
            try:
                rp = os.path.realpath(p)
            except OSError:
                rp = p
            if rp in seen or not os.path.isfile(p):
                continue
            seen.add(rp)
            out.append(p)
    return sorted(out)


def collect_dark_genes_from_analysis_dir(
    analysis_dir: str,
    sample_name: str = "",
) -> Dict[str, Any]:
    """
    Scan analysis_dir for dark-gene pipeline summary artifacts.

    Returns a dict suitable for result.json[\"dark_genes\"] and report.json.
    """
    base = (analysis_dir or "").strip()
    if not base or not os.path.isdir(base):
        return {
            "status": "not_found",
            "message": "analysis_dir missing or not a directory",
        }

    summary_paths = _first_existing_glob(base, _SUMMARY_GLOBS)
    detailed_paths = _first_existing_glob(base, _DETAILED_GLOBS)

    if not summary_paths and not detailed_paths:
        return {
            "status": "not_found",
            "message": "No dark_genes summary/detailed report under analysis output",
            "searched_root": os.path.abspath(base),
        }

    summary_text = ""
    if summary_paths:
        summary_text = _read_text_safe(summary_paths[0], _MAX_READ_BYTES)

    detailed_excerpt = ""
    if detailed_paths:
        detailed_excerpt = _read_text_safe(detailed_paths[0], min(120_000, _MAX_READ_BYTES))

    status = "found" if summary_text else "partial"

    return {
        "status": status,
        "sample_name": (sample_name or "").strip(),
        "summary_file": os.path.basename(summary_paths[0]) if summary_paths else None,
        "detailed_file": os.path.basename(detailed_paths[0]) if detailed_paths else None,
        "summary_paths": [os.path.relpath(p, base) for p in summary_paths],
        "detailed_paths": [os.path.relpath(p, base) for p in detailed_paths],
        "summary_text": summary_text,
        "detailed_text": detailed_excerpt,
    }


def dark_genes_for_pdf(report_block: Dict[str, Any]) -> Dict[str, Any]:
    """
    Subset + truncate for customer PDF (plain text; Jinja escapes HTML).
    """
    if not report_block or report_block.get("status") == "not_found":
        return {}
    if report_block.get("status") == "error":
        msg = (report_block.get("message") or "unknown error")[:4000]
        return {
            "status": "error",
            "report_summary": f"Supplementary dark-gene analysis metadata could not be loaded: {msg}",
        }
    summary = (report_block.get("summary_text") or "").strip()
    if len(summary) > 12_000:
        summary = summary[:12_000] + "\n\n[… truncated for PDF …]"
    detailed = (report_block.get("detailed_text") or "").strip()
    if len(detailed) > 8_000:
        detailed = detailed[:8_000] + "\n\n[… truncated …]"
    return {
        "status": report_block.get("status"),
        "report_summary": summary,
        "report_detailed": detailed if detailed else None,
    }
