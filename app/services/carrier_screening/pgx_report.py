"""
PGx (PharmCAT) pipeline outputs under <analysis_or_output>/pgx/ → result.json + customer PDF.

Expected layout (gx-exome):
  pgx/pgx_meta.json
  pgx/pgx_summary.txt
  pgx/pgx_result.json          (artifact pointers; may be large — not embedded wholesale)
  pgx/<sample>_pgx.report.html (full PharmCAT HTML — not embedded in PDF; use summary + meta)

The customer PDF uses the plain-text summary (human-readable, WeasyPrint-safe).

Portal also gets ``gene_results``: one row per gene from ``pgx_result.json`` → ``phenotype``
(PharmCAT ``geneReports``), with ``reviewer_confirmed`` / ``reviewer_comment`` merged on reprocess.
"""

from __future__ import annotations

import glob
import html
import json
import logging
import os
from typing import AbstractSet, Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# Mirrors gx-exome `pgx.nf` PYSUMMARY (actionable vs normal)
_RISK_FUNCTIONS = frozenset({"no function", "decreased function", "unfavorable response allele"})
_SKIP_PHENOTYPES = frozenset({"no result", "n/a", ""})


def extract_gene_results_from_phenotype(phenotype: Any) -> List[Dict[str, Any]]:
    """
    Build one row per gene from PharmCAT phenotype JSON (``geneReports`` / recommendation diplotypes).

    Same selection rules as gx-exome ``pgx.nf`` summary script so the portal table matches the text summary.
    """
    if not isinstance(phenotype, dict):
        return []
    rows: List[Dict[str, Any]] = []
    seen_genes: set = set()
    for source in ("CPIC", "DPWG"):
        gene_reports = phenotype.get("geneReports")
        if not isinstance(gene_reports, dict):
            continue
        reports = gene_reports.get(source)
        if not isinstance(reports, dict):
            continue
        for gene_name in sorted(reports.keys()):
            if gene_name in seen_genes:
                continue
            gene_data = reports.get(gene_name)
            if not isinstance(gene_data, dict):
                continue
            call_src = (gene_data.get("callSource") or "").strip()
            rec_dips = gene_data.get("recommendationDiplotypes")
            if not isinstance(rec_dips, list):
                continue
            for dip in rec_dips:
                if not isinstance(dip, dict):
                    continue
                a1 = dip.get("allele1") if isinstance(dip.get("allele1"), dict) else {}
                a2 = dip.get("allele2") if isinstance(dip.get("allele2"), dict) else {}
                n1 = (a1.get("name") or "").strip()
                n2 = (a2.get("name") or "").strip()
                fn1 = (a1.get("function") or "").strip()
                fn2 = (a2.get("function") or "").strip()
                phenotypes = dip.get("phenotypes")
                if not isinstance(phenotypes, list):
                    phenotypes = []
                activity = dip.get("activityScore")
                diplotype = f"{n1}/{n2}" if n2 else n1
                phenotype_str = ", ".join(p for p in phenotypes if p) if phenotypes else ""
                low = phenotype_str.lower()
                if low in _SKIP_PHENOTYPES:
                    continue
                functions = [f for f in (fn1, fn2) if f]
                has_risk = any((f or "").lower() in _RISK_FUNCTIONS for f in functions)
                cat = "actionable" if has_risk else "normal"
                rows.append(
                    {
                        "gene": gene_name,
                        "guideline_source": source,
                        "diplotype": diplotype,
                        "phenotype": phenotype_str,
                        "activity_score": activity,
                        "allele1_function": fn1,
                        "allele2_function": fn2,
                        "call_source": call_src,
                        "category": cat,
                    }
                )
                seen_genes.add(gene_name)
                break
    rows.sort(
        key=lambda r: (
            0 if r.get("category") == "actionable" else 1,
            str(r.get("gene") or ""),
        )
    )
    return rows


def filter_pgx_gene_results_by_panel(
    gene_results: List[Dict[str, Any]],
    panel_genes: AbstractSet[str],
) -> List[Dict[str, Any]]:
    """
    Keep only PGx rows whose ``gene`` is in the order interpretation set (WES panel + extras).

    Same gene universe as variant post-filter: ``interpretation_gene_set_for_job`` / ``panel_interpretation_genes``.
    When ``panel_genes`` is empty, returns ``gene_results`` unchanged (caller should not pass an empty set to mean “filter”).
    """
    if not panel_genes:
        return list(gene_results)
    up = {str(g).strip().upper() for g in panel_genes if g is not None and str(g).strip()}
    if not up:
        return list(gene_results)
    out: List[Dict[str, Any]] = []
    for row in gene_results:
        if not isinstance(row, dict):
            continue
        g = row.get("gene")
        if g is None or str(g).strip().upper() not in up:
            continue
        out.append(row)
    return out


def merge_pgx_gene_reviews(
    fresh_rows: List[Dict[str, Any]],
    previous_pgx: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """After reprocess, keep ``reviewer_confirmed`` / ``reviewer_comment`` per ``gene``."""
    prev_by_gene: Dict[str, Dict[str, Any]] = {}
    if isinstance(previous_pgx, dict):
        for row in previous_pgx.get("gene_results") or []:
            if isinstance(row, dict) and row.get("gene"):
                prev_by_gene[str(row["gene"])] = row
    out: List[Dict[str, Any]] = []
    for row in fresh_rows:
        if not isinstance(row, dict):
            continue
        g = row.get("gene")
        merged = dict(row)
        pr = prev_by_gene.get(str(g)) if g is not None else None
        if pr:
            merged["reviewer_confirmed"] = bool(pr.get("reviewer_confirmed", False))
            merged["reviewer_comment"] = (pr.get("reviewer_comment") or "")[:4000]
        else:
            merged.setdefault("reviewer_confirmed", False)
            merged.setdefault("reviewer_comment", "")
        out.append(merged)
    return out


def merge_pgx_portal_review(
    fresh: Dict[str, Any],
    previous: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    When ``generate_result_json`` rewrites ``pgx``, keep ``portal_review`` from the prior
    ``result.json`` (reviewer notes / reviewed flag from the portal).
    """
    out = dict(fresh) if isinstance(fresh, dict) else {}
    if not isinstance(previous, dict):
        return out
    pr = previous.get("portal_review")
    if isinstance(pr, dict) and pr:
        out["portal_review"] = dict(pr)
    return out


def collect_pgx_from_analysis_dir(root: str, sample_name: str) -> Dict[str, Any]:
    """
    Scan ``root/pgx/`` for PharmCAT meta + summary text.

    Returns a dict suitable for ``result.json["pgx"]``.
    """
    root = (root or "").strip()
    if not root or not os.path.isdir(root):
        return {
            "status": "not_found",
            "message": "analysis/output root is missing or not a directory",
            "searched": root,
        }

    pgx_dir = os.path.join(root, "pgx")
    if not os.path.isdir(pgx_dir):
        return {
            "status": "not_found",
            "message": f"No pgx directory under {root}",
            "searched_root": os.path.abspath(root),
        }

    meta: Dict[str, Any] = {}
    meta_path = os.path.join(pgx_dir, "pgx_meta.json")
    if os.path.isfile(meta_path):
        try:
            with open(meta_path, "r", encoding="utf-8") as f:
                meta = json.load(f)
        except Exception as e:
            logger.warning("[pgx] could not read %s: %s", meta_path, e)

    summary_text = ""
    summary_path = os.path.join(pgx_dir, "pgx_summary.txt")
    if os.path.isfile(summary_path):
        try:
            with open(summary_path, "r", encoding="utf-8") as f:
                summary_text = f.read()
        except Exception as e:
            logger.warning("[pgx] could not read %s: %s", summary_path, e)

    reporter_html = _first_basename(pgx_dir, "*_pgx.report.html")
    result_json_name: Optional[str] = None
    rp = os.path.join(pgx_dir, "pgx_result.json")
    if os.path.isfile(rp):
        result_json_name = "pgx_result.json"

    exit_ok = str((meta or {}).get("exit_status") or "").lower() in ("success", "ok", "0")
    if not summary_text.strip():
        if exit_ok and meta:
            summary_text = (
                f"PGx completed ({(meta.get('tool_version') or 'PharmCAT')}). "
                "Summary file was empty; see full PharmCAT HTML in pipeline pgx/ output."
            )
        else:
            return {
                "status": "error",
                "message": "pgx_summary.txt missing or empty",
                "meta": meta,
                "pgx_dir": os.path.abspath(pgx_dir),
            }

    out: Dict[str, Any] = {
        "status": "ok",
        "summary_text": summary_text,
        "meta": meta,
        "artifacts": {
            "pgx_summary_txt": "pgx/pgx_summary.txt",
            "pgx_meta_json": "pgx/pgx_meta.json",
            "reporter_html_basename": reporter_html,
            "pgx_result_json": result_json_name,
        },
        "pgx_dir": os.path.abspath(pgx_dir),
        "gene_results": [],
    }

    rp_full = os.path.join(pgx_dir, "pgx_result.json")
    if os.path.isfile(rp_full):
        try:
            with open(rp_full, "r", encoding="utf-8") as f:
                bundle = json.load(f)
            phen = bundle.get("phenotype") if isinstance(bundle, dict) else None
            if isinstance(phen, dict):
                out["gene_results"] = extract_gene_results_from_phenotype(phen)
        except Exception as e:
            logger.warning("[pgx] could not parse gene rows from %s: %s", rp_full, e)

    return out


def _first_basename(dir_path: str, pattern: str) -> Optional[str]:
    paths = sorted(glob.glob(os.path.join(dir_path, pattern)))
    if not paths:
        return None
    return os.path.basename(paths[0])


def pgx_for_pdf(pgx: Dict[str, Any]) -> Dict[str, Any]:
    """
    Subset + HTML fragment for WeasyPrint (written into report.json used by Jinja).
    """
    if not isinstance(pgx, dict):
        return {}
    st = pgx.get("status")
    if st == "not_found":
        return {}

    meta = pgx.get("meta") if isinstance(pgx.get("meta"), dict) else {}
    tool_v = (meta.get("tool_version") or meta.get("tool") or "PharmCAT").strip()

    if st == "error":
        msg = html.escape(str(pgx.get("message") or "PGx error"))
        err_out: Dict[str, Any] = {
            **pgx,
            "summary_for_pdf_html": (
                f'<p class="detail-text" style="color:#b91c1c;font-size:9pt;">{msg}</p>'
            ),
            "tool_version_line": tool_v,
        }
        err_out.pop("portal_review", None)
        err_out.pop("gene_results", None)
        return err_out

    body = pgx.get("summary_text") or ""
    esc = html.escape(body)
    foot = (
        f'<p class="muted" style="font-size:7.5pt;margin-top:10px;line-height:1.4;">'
        f"Source: PharmCAT summary text. Full interactive report: "
        f"<code>pgx/</code> output (<code>*_pgx.report.html</code>). {html.escape(tool_v)}"
        f"</p>"
    )
    block = (
        f'<pre class="pgx-summary" style="white-space:pre-wrap;font-family:ui-monospace,monospace;'
        f"font-size:8.5pt;line-height:1.35;border:1px solid #cbd5e1;border-radius:8px;"
        f'padding:12px;background:#f8fafc;">{esc}</pre>{foot}'
    )
    out = dict(pgx)
    out["summary_for_pdf_html"] = block
    out["tool_version_line"] = tool_v
    out.pop("portal_review", None)
    out.pop("gene_results", None)
    return out


def sanitize_pgx_payload_for_pdf_render(report_data: Dict[str, Any]) -> None:
    """Ensure ``report_data['pgx']`` has ``summary_for_pdf_html`` when raw summary exists."""
    pgx = report_data.get("pgx")
    if not isinstance(pgx, dict):
        return
    if (pgx.get("summary_for_pdf_html") or "").strip():
        return
    merged = pgx_for_pdf(pgx)
    if merged:
        combined = {**pgx, **merged}
        combined.pop("portal_review", None)
        combined.pop("gene_results", None)
        report_data["pgx"] = combined
