"""
Gene ↔ disease association helpers.

1. Panel JSON via DiseaseGeneMapper (primary in VariantAnnotator during VCF annotation).
2. SQLite / Gemini for disorder names: :mod:`gene_knowledge_db` — **only** when generating a report for **confirmed** variants (see ``enrich_confirmed_variants_for_report``), or via CLI ``populate_gene_knowledge_db.py``.
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from .gene_knowledge_db import (
    diseases_from_gene_knowledge_sqlite,
    enrich_confirmed_variants_for_report,
    fetch_gene_knowledge_via_gemini,
    get_or_fetch_gene_diseases,
    init_gene_knowledge_database,
)

logger = logging.getLogger(__name__)

__all__ = [
    "diseases_from_gene_knowledge_sqlite",
    "enrich_confirmed_variants_for_report",
    "merge_diseases",
    "fetch_remote_gemini_sample",
    "get_or_fetch_gene_diseases",
    "init_gene_knowledge_database",
]


def merge_diseases(panel: List[Dict[str, Any]], fallback: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    if panel:
        return panel
    return fallback


def fetch_remote_gemini_sample(
    gene: str,
    api_key: str,
    model: str = "gemini-2.5-flash",
) -> Optional[Dict[str, Any]]:
    """
    One-off Gemini lookup (no DB). Returns annotator-shaped disease dict or None.
    """
    flat, _err = fetch_gene_knowledge_via_gemini(gene, api_key, model=model)
    if not flat:
        return None
    dn = (flat.get("disorder") or "").strip()
    if not dn:
        return None
    omim = (flat.get("omim_number") or "").strip().replace("OMIM:", "")
    return {
        "name": dn,
        "disease_name": dn,
        "omim_id": omim,
        "inheritance": (flat.get("inheritance") or "").strip(),
        "source": "gemini_search",
    }
