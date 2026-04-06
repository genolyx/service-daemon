"""
Portal: gene vs ordered WES / carrier panel — target BED intervals and optional per-gene depth.

Uses ``job.params`` + panel catalog (``wes_panel_id``) to resolve ``disease_bed``, then lists BED
rows whose 4th column matches the HGNC symbol. Optionally scans pipeline output for small
``*gene*coverage*`` / ``*coverage*by*gene*`` TSV-style sidecars (Twist / custom QC).
"""

from __future__ import annotations

import csv
import glob
import gzip
import json
import logging
import os
import re
from typing import Any, Dict, List, Optional, Tuple

from ..models import Job
from .wes_panels import (
    _resolve_repo_path,
    get_panel_by_id,
    interpretation_gene_set_for_job,
)

logger = logging.getLogger(__name__)

_CARRIER_LIKE = frozenset({"carrier_screening", "whole_exome", "health_screening"})

_GENE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9-]{0,24}$")

# mosdepth per-base: chrom, start, end, depth (BED intervals, depth constant on each segment)
_MOSDEPTH_PER_BASE_PATTERNS = (
    "*mosdepth*.per-base.bed.gz",
    "*mosdepth*per-base*.bed.gz",
    "*.mosdepth.per-base.bed.gz",
)


def normalize_gene_symbol(raw: str) -> str:
    return (raw or "").strip().upper()


def _artifact_path_ok(path: str) -> bool:
    norm = path.replace("\\", "/").lower()
    return "/env/" not in norm and "/viz_env/" not in norm


def _chrom_name_variants(chrom: str) -> List[str]:
    """Tabix may index chr1 vs 1 — try both."""
    c = (chrom or "").strip()
    if not c:
        return []
    out: List[str] = [c]
    if c.lower().startswith("chr"):
        out.append(c[3:])
    else:
        out.append("chr" + c)
    seen: set = set()
    uniq: List[str] = []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq


def overlap_bp(g0: int, g1: int, seg0: int, seg1: int) -> int:
    """Overlap length between [g0,g1) and [seg0,seg1) (half-open)."""
    a = max(seg0, g0)
    b = min(seg1, g1)
    return max(0, b - a)


def _find_mosdepth_per_base_bed(roots: List[str]) -> Optional[str]:
    """Prefer indexed mosdepth ``*.per-base.bed.gz`` next to ``*.tbi``."""
    candidates: List[Tuple[float, str]] = []
    for root in roots:
        if not root or not os.path.isdir(root):
            continue
        for pat in _MOSDEPTH_PER_BASE_PATTERNS:
            try:
                for fp in glob.glob(os.path.join(root, "**", pat), recursive=True):
                    if not os.path.isfile(fp) or not _artifact_path_ok(fp):
                        continue
                    if not fp.endswith(".gz"):
                        continue
                    if not os.path.isfile(fp + ".tbi"):
                        continue
                    try:
                        candidates.append((os.path.getmtime(fp), fp))
                    except OSError:
                        continue
            except Exception:
                continue
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def _pct_bases_ge_depths_from_per_base(
    regions: List[Dict[str, Any]],
    per_base_gz: str,
    thresholds: Tuple[int, ...] = (10, 20),
) -> Optional[Dict[str, Any]]:
    """
    Intersect panel BED intervals for the gene with mosdepth per-base segments; compute % of
    assessed bases at or above each depth threshold.
    """
    if not regions or not per_base_gz or not os.path.isfile(per_base_gz + ".tbi"):
        return None
    try:
        import pysam
    except ImportError:
        logger.warning("[gene_panel_coverage] pysam not available for per-base depth")
        return None

    tb = None
    try:
        tb = pysam.TabixFile(per_base_gz)
    except Exception as e:
        logger.debug("[gene_panel_coverage] TabixFile open failed %s: %s", per_base_gz, e)
        return None

    # total bases per threshold (same total for all)
    total_bp = 0
    ge_counts = {t: 0 for t in thresholds}

    try:
        for reg in regions:
            g0 = int(reg.get("start", 0))
            g1 = int(reg.get("end", 0))
            if g1 <= g0:
                continue
            chrom = str(reg.get("chrom") or "").strip()
            fetched = False
            for ctry in _chrom_name_variants(chrom):
                try:
                    for row in tb.fetch(ctry, g0, g1):
                        parts = row.strip().split("\t")
                        if len(parts) < 4:
                            continue
                        try:
                            s = int(parts[1])
                            e = int(parts[2])
                            d = float(parts[3])
                        except (ValueError, TypeError):
                            continue
                        L = overlap_bp(g0, g1, s, e)
                        if L <= 0:
                            continue
                        total_bp += L
                        for t in thresholds:
                            if d >= t:
                                ge_counts[t] += L
                    fetched = True
                    break
                except Exception as e:
                    logger.debug("[gene_panel_coverage] tabix fetch %s %s-%s: %s", ctry, g0, g1, e)
                    continue
            if not fetched:
                continue
    finally:
        try:
            if tb is not None:
                tb.close()
        except Exception:
            pass

    if total_bp <= 0:
        return None

    out: Dict[str, Any] = {
        "total_bases_assessed": int(total_bp),
        "source_file": per_base_gz,
        "method": "mosdepth_per_base_tabix",
    }
    for t in thresholds:
        pct = round(100.0 * float(ge_counts[t]) / float(total_bp), 2)
        out[f"pct_bases_ge_{t}x"] = pct
    return out


def _wes_panel_id_from_job(job: Job) -> str:
    p = job.params or {}
    pid = p.get("wes_panel_id")
    if not pid and isinstance(p.get("carrier"), dict):
        pid = p["carrier"].get("wes_panel_id")
    return (str(pid).strip() if pid else "") or ""


def _explicit_disease_bed_path(job: Job, panel: Optional[Dict[str, Any]]) -> Tuple[Optional[str], str]:
    """``disease_bed`` from job params or panel catalog only (no backbone)."""
    p = job.params or {}
    db = (p.get("disease_bed") or "").strip()
    if db:
        db_abs = os.path.normpath(db)
        if os.path.isfile(db_abs):
            return db_abs, "job.params.disease_bed"
    if panel:
        rel = str(panel.get("disease_bed") or "").strip()
        if rel:
            abs_p = _resolve_repo_path(rel)
            if os.path.isfile(abs_p):
                return abs_p, "panel_catalog.disease_bed"
    return None, ""


def _effective_backbone_bed_path(job: Job, panel: Optional[Dict[str, Any]]) -> Tuple[Optional[str], str]:
    """Capture / exome target BED (HGNC in col4 for Twist-style kits)."""
    p = job.params or {}
    bb = (p.get("backbone_bed") or "").strip()
    if bb:
        bb_abs = os.path.normpath(bb)
        if os.path.isfile(bb_abs):
            return bb_abs, "job.params.backbone_bed"
    if panel:
        rel = str(panel.get("backbone_bed") or "").strip()
        if rel:
            abs_p = _resolve_repo_path(rel)
            if os.path.isfile(abs_p):
                return abs_p, "panel_catalog.backbone_bed"
    try:
        from ..config import settings

        for opt, label in (
            (getattr(settings, "carrier_default_backbone_bed", None), "settings.carrier_default_backbone_bed"),
            (getattr(settings, "wes_panel_gene_source_bed", None), "settings.wes_panel_gene_source_bed"),
        ):
            if not opt or not str(opt).strip():
                continue
            ap = os.path.normpath(str(opt).strip())
            if os.path.isfile(ap):
                return ap, label
    except Exception as e:
        logger.debug("[gene_panel_coverage] backbone settings fallback: %s", e)
    return None, ""


def _effective_disease_bed_path(job: Job, panel: Optional[Dict[str, Any]]) -> Tuple[Optional[str], str]:
    """
    Last-resort: any BED path for gene-interval lookup when explicit disease + backbone yield no rows.

    Order: explicit disease → panel disease → backbone chain → disease defaults → gene_bed.
    """
    dis_p, dis_s = _explicit_disease_bed_path(job, panel)
    if dis_p:
        return dis_p, dis_s
    bb_p, bb_s = _effective_backbone_bed_path(job, panel)
    if bb_p:
        return bb_p, bb_s
    try:
        from ..config import settings

        for opt, label in (
            (getattr(settings, "carrier_default_disease_bed", None), "settings.carrier_default_disease_bed"),
            (getattr(settings, "gene_bed", None), "settings.gene_bed"),
        ):
            if not opt or not str(opt).strip():
                continue
            ap = os.path.normpath(str(opt).strip())
            if os.path.isfile(ap):
                return ap, label
    except Exception as e:
        logger.debug("[gene_panel_coverage] settings BED fallback: %s", e)
    return None, ""


def _bed_col4_matches_gene(name_field: str, gene: str) -> bool:
    """
    True if *gene* (HGNC symbol) is named in BED column 4.

    Plain BEDs use ``PAH`` alone. Twist / vendor BEDs often use composite labels such as
    ``PAH;NM_000277.3;ENST00000307000.7;ClinID-…`` — require token match, not full-string equality.
    """
    g = normalize_gene_symbol(gene)
    if not g:
        return False
    raw = (name_field or "").strip()
    if not raw:
        return False
    if raw.upper() == g:
        return True
    for sep in (";", ",", "|"):
        if sep not in raw:
            continue
        for tok in raw.split(sep):
            t = tok.strip()
            if not t:
                continue
            if t.upper() == g:
                return True
    return False


def _bed_intervals_for_gene(bed_path: str, gene: str) -> List[Dict[str, Any]]:
    g = gene.upper()
    rows: List[Dict[str, Any]] = []
    if not bed_path or not os.path.isfile(bed_path):
        return rows
    opener = gzip.open if bed_path.endswith(".gz") else open
    try:
        with opener(bed_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                name = (parts[3].strip() if len(parts) >= 4 else "") or ""
                if not _bed_col4_matches_gene(name, g):
                    continue
                chrom = parts[0].strip()
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except (ValueError, TypeError):
                    continue
                bp = max(0, end - start)
                rows.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "name": name,
                        "length_bp": bp,
                    }
                )
    except OSError as e:
        logger.warning("[gene_panel_coverage] Cannot read BED %s: %s", bed_path, e)
    return rows


def _coverage_search_roots(job: Job) -> List[str]:
    roots: List[str] = []
    seen: set = set()

    def add(p: Optional[str]) -> None:
        if not p or not str(p).strip():
            return
        path = str(p).strip()
        if not os.path.isdir(path):
            return
        try:
            r = os.path.realpath(path)
        except OSError:
            return
        if r in seen:
            return
        seen.add(r)
        roots.append(os.path.abspath(path))

    if job.service_code in _CARRIER_LIKE:
        try:
            from . import get_plugin

            pl = get_plugin(job.service_code)
            if pl is not None and hasattr(pl, "dark_genes_search_roots"):
                for x in pl.dark_genes_search_roots(job):
                    add(x)
        except Exception as e:
            logger.debug("[gene_panel_coverage] dark_genes_search_roots: %s", e)

    add(getattr(job, "analysis_dir", None))
    add(getattr(job, "output_dir", None))
    if job.service_code in _CARRIER_LIKE:
        try:
            from .carrier_screening.plugin import carrier_report_output_dir

            add(carrier_report_output_dir(job))
        except Exception:
            pass
    return roots


def _parse_pct_cell(raw: str) -> Optional[float]:
    """Interpret table cell as percent (accepts 0–1 fraction or 0–100)."""
    if raw is None or not str(raw).strip():
        return None
    try:
        v = float(str(raw).strip().replace(",", "").replace("%", ""))
    except ValueError:
        return None
    if 0 <= v <= 1.0:
        return round(v * 100.0, 2)
    return round(v, 2)


def _hdr_find_pct_col(hdr: List[str], subs: Tuple[str, ...]) -> Optional[int]:
    for i, h in enumerate(hdr):
        h0 = h.replace(" ", "_").replace("-", "_").lower()
        for s in subs:
            if s in h0:
                return i
    return None


def _read_carrier_qc_summary(job: Job) -> Dict[str, Any]:
    if job.service_code not in _CARRIER_LIKE:
        return {}
    try:
        from .carrier_screening.plugin import carrier_result_json_path

        path = carrier_result_json_path(job)
    except Exception:
        return {}
    if not path or not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        qc = data.get("qc_summary")
        return qc if isinstance(qc, dict) else {}
    except Exception as e:
        logger.debug("[gene_panel_coverage] qc_summary read failed: %s", e)
        return {}


def _parse_gene_depth_file(path: str, gene: str) -> Optional[Dict[str, Any]]:
    """Best-effort: TSV/CSV with a gene column and a numeric depth / mean column."""
    g = gene.upper()
    try:
        if os.path.getsize(path) > 800 * 1024:
            return None
    except OSError:
        return None
    try:
        with open(path, "r", errors="replace", newline="") as f:
            sample = f.read(64 * 1024)
            f.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            except Exception:
                dialect = csv.excel_tab
            reader = csv.reader(f, dialect)
            rows = list(reader)
    except Exception:
        try:
            with open(path, "r", errors="replace") as f:
                rows = [ln.rstrip("\n").split("\t") for ln in f.readlines()]
        except OSError:
            return None

    if not rows:
        return None

    header = [str(c).strip().lower() for c in rows[0]]
    gene_keys = ("gene", "gene_symbol", "symbol", "hgnc", "hgnc_id")
    depth_keys = (
        "mean_coverage",
        "mean_depth",
        "avg_depth",
        "meandepth",
        "depth",
        "coverage",
        "cov",
        "avg_coverage",
    )

    def col_idx(keys: Tuple[str, ...], hdr: List[str]) -> Optional[int]:
        for i, h in enumerate(hdr):
            h0 = h.replace(" ", "_").replace("-", "_")
            for k in keys:
                if k == h0 or (k in h0 and len(h0) < 40):
                    return i
        return None

    gi = col_idx(gene_keys, header)
    di = col_idx(depth_keys, header)
    pi10 = _hdr_find_pct_col(header, ("pct_bases_10x", "pct_10x", "bases_ge_10", "ge_10x", ">=_10x"))
    pi20 = _hdr_find_pct_col(header, ("pct_bases_20x", "pct_20x", "bases_ge_20", "ge_20x", ">=_20x"))

    data_rows = rows[1:] if gi is not None and di is not None else rows

    for parts in data_rows:
        if not parts:
            continue
        if gi is not None and di is not None:
            idx_needed = [gi, di, pi10, pi20]
            mx = max(i for i in idx_needed if i is not None)
            if len(parts) <= mx:
                continue
            sym = str(parts[gi]).strip().upper()
            if sym != g:
                continue
            raw = str(parts[di]).strip().replace(",", "")
            try:
                val = float(raw)
            except ValueError:
                continue
            out: Dict[str, Any] = {"mean_coverage": val, "method": "tabular_file"}
            if pi10 is not None:
                p = _parse_pct_cell(str(parts[pi10]))
                if p is not None:
                    out["pct_bases_ge_10x"] = p
            if pi20 is not None:
                p = _parse_pct_cell(str(parts[pi20]))
                if p is not None:
                    out["pct_bases_ge_20x"] = p
            return out
        else:
            # Two-column: GENE<TAB>depth
            if len(parts) >= 2 and str(parts[0]).strip().upper() == g:
                raw = str(parts[1]).strip().replace(",", "")
                try:
                    val = float(raw)
                    return {"mean_coverage": val, "method": "two_column"}
                except ValueError:
                    continue
    return None


def _scan_per_gene_depth(roots: List[str], gene: str) -> Optional[Dict[str, Any]]:
    patterns = (
        "*gene*coverage*.txt",
        "*gene*coverage*.tsv",
        "*coverage*by*gene*.txt",
        "*coverage*by*gene*.tsv",
        "*per*gene*coverage*.txt",
        "*per*gene*coverage*.tsv",
    )
    candidates: List[Tuple[float, str]] = []
    for root in roots:
        if not root or not os.path.isdir(root):
            continue
        for pat in patterns:
            try:
                for fp in glob.glob(os.path.join(root, "**", pat), recursive=True):
                    if not os.path.isfile(fp) or not _artifact_path_ok(fp):
                        continue
                    try:
                        candidates.append((os.path.getmtime(fp), fp))
                    except OSError:
                        continue
            except Exception:
                continue
    candidates.sort(key=lambda x: x[0], reverse=True)

    seen: set = set()
    for _, fp in candidates:
        if fp in seen:
            continue
        seen.add(fp)
        if len(seen) > 48:
            break
        hit = _parse_gene_depth_file(fp, gene)
        if hit:
            out = dict(hit)
            out["source_file"] = fp
            return out
    return None


def build_gene_panel_coverage_report(job: Job, gene_raw: str) -> Dict[str, Any]:
    """
    Build a JSON-serializable report for the portal.
    """
    gene = normalize_gene_symbol(gene_raw)
    if not gene or not _GENE_RE.match(gene_raw.strip() or ""):
        raise ValueError("Invalid gene symbol")

    wpid = _wes_panel_id_from_job(job)
    panel = get_panel_by_id(wpid) if wpid else None
    panel_label = (panel.get("label") or wpid) if panel else None

    interp = interpretation_gene_set_for_job(job)
    in_set = gene in interp

    dis_path, dis_src = _explicit_disease_bed_path(job, panel)
    bb_path, bb_src = _effective_backbone_bed_path(job, panel)
    r_bb = _bed_intervals_for_gene(bb_path or "", gene) if bb_path else []
    r_dis = _bed_intervals_for_gene(dis_path or "", gene) if dis_path else []

    clinical_disease_bed_path: Optional[str] = None
    clinical_disease_bed_source: Optional[str] = None

    if r_bb:
        regions = r_bb
        bed_path, bed_source = bb_path, bb_src
        if dis_path and dis_path != bb_path:
            clinical_disease_bed_path, clinical_disease_bed_source = dis_path, dis_src
    elif r_dis:
        regions = r_dis
        bed_path, bed_source = dis_path, dis_src
    else:
        bed_path, bed_source = _effective_disease_bed_path(job, panel)
        regions = _bed_intervals_for_gene(bed_path or "", gene) if bed_path else []

    total_bp = sum(int(r.get("length_bp") or 0) for r in regions)

    qc = _read_carrier_qc_summary(job)
    cov = qc.get("coverage") if isinstance(qc.get("coverage"), dict) else {}
    global_cov: Dict[str, Any] = {}
    if isinstance(cov, dict):
        for k in ("mean_coverage", "min_coverage", "max_coverage"):
            if cov.get(k) is not None:
                global_cov[k] = cov[k]
        for k, v in cov.items():
            if isinstance(k, str) and k.startswith("pct_bases_") and v is not None:
                global_cov[k] = v

    roots = _coverage_search_roots(job)
    per_gene = _scan_per_gene_depth(roots, gene)

    gene_depth_thresholds: Optional[Dict[str, Any]] = None
    if regions:
        pb = _find_mosdepth_per_base_bed(roots)
        if pb:
            gene_depth_thresholds = _pct_bases_ge_depths_from_per_base(regions, pb)

    if gene_depth_thresholds is None and per_gene:
        p10 = per_gene.get("pct_bases_ge_10x")
        p20 = per_gene.get("pct_bases_ge_20x")
        if p10 is not None or p20 is not None:
            gene_depth_thresholds = {
                "method": "gene_qc_sidecar",
                "source_file": per_gene.get("source_file"),
            }
            if p10 is not None:
                gene_depth_thresholds["pct_bases_ge_10x"] = p10
            if p20 is not None:
                gene_depth_thresholds["pct_bases_ge_20x"] = p20

    notes: List[str] = []
    if bed_path and regions:
        if clinical_disease_bed_path:
            notes.append(
                "Depth percentages use **capture (backbone) target** intervals for this gene (bases typically sequenced at high depth on this exome). "
                "A separate clinical/disease panel BED may span more bases (e.g. intronic or off-capture), where depth is often low — that mismatch can make disease-BED-only stats look worse than run-wide exome QC."
            )
        else:
            src = (bed_source or "").lower()
            if "backbone" in src:
                notes.append(
                    "Regions are from the capture/backbone BED (column 4 matches this gene); "
                    "no separate disease BED was set on this order."
                )
            else:
                notes.append(
                    "Regions are rows from the panel disease BED whose 4th column matches this gene (design / capture target)."
                )
    elif bed_path and not regions:
        notes.append(
            "A BED path resolved, but no rows use this gene name in column 4 — the gene may be listed only in interpretation_genes, or labels differ from HGNC symbols."
        )
    if not bed_path:
        notes.append("No BED resolved for this order/panel — interval list unavailable.")

    if gene_depth_thresholds and gene_depth_thresholds.get("method") == "mosdepth_per_base_tabix":
        notes.append(
            "% of bases ≥10× / ≥20× is from mosdepth per-base depth segments intersected with the BED intervals above (half-open coordinates)."
        )
    elif regions and not gene_depth_thresholds:
        notes.append(
            "Gene-level % bases at 10×/20× not computed: publish an indexed mosdepth per-base file "
            "(*mosdepth*.per-base.bed.gz + .tbi) under analysis/output, or a sidecar with pct_bases_10x / pct_bases_20x."
        )

    return {
        "gene": gene,
        "order_id": job.order_id,
        "service_code": job.service_code,
        "wes_panel_id": wpid or None,
        "panel_label": panel_label,
        "in_interpretation_set": in_set,
        "interpretation_gene_count": len(interp),
        "disease_bed_path": bed_path,
        "disease_bed_source": bed_source or None,
        "clinical_disease_bed_path": clinical_disease_bed_path,
        "clinical_disease_bed_source": clinical_disease_bed_source,
        "bed_regions": regions,
        "total_target_bp": total_bp,
        "panel_bed_region_count": len(regions),
        "global_qc_coverage": global_cov,
        "per_gene_depth": per_gene,
        "gene_depth_thresholds": gene_depth_thresholds,
        "notes": notes,
    }
