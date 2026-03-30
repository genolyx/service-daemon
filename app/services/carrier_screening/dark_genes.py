"""
Collect dark-gene pipeline outputs (unified Nextflow outdir = analysis_dir).

Looks for summary reports under summary/ or results/summary/ (dark_gene_pipeline layout),
or *_summary_report.txt / *_detailed_report.txt in the analysis outdir root (some runs publish flat).
"""

from __future__ import annotations

import glob
import html
import logging
import os
import re
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def _coerce_approved_bool(val: Any) -> bool:
    """
    Normalize ``approved`` from JSON/API. Python ``bool("false")`` is True — strings must
    be interpreted explicitly.
    """
    if val is True:
        return True
    if val is False or val is None:
        return False
    if isinstance(val, (int, float)):
        return bool(val)
    if isinstance(val, str):
        s = val.strip().lower()
        if s in ("", "0", "false", "no", "n", "off"):
            return False
        if s in ("1", "true", "yes", "y", "on"):
            return True
        return False
    return False


def _coerce_risk_level(val: Any) -> str:
    """``low`` (green accent) vs ``high`` (red accent); explicit strings only."""
    s = str(val).strip().lower() if val is not None else ""
    if s == "low":
        return "low"
    return "high"


def _infer_pipeline_section_high_risk(sec: Dict[str, Any]) -> bool:
    """
    True when the pipeline marked this block as risky (``WARNING:`` line or warning-kind section).
    Used when ``section_reviews[i].risk`` is absent — default is **low** unless this is true.
    """
    if not sec or not isinstance(sec, dict):
        return False
    k = _infer_section_kind(sec)
    if k == "warning":
        return True
    tl = (sec.get("title") or "").strip().upper()
    if "WARNING" in tl and "QUALITY" not in tl:
        return True
    body = (sec.get("body") or "")
    if re.search(r"(?im)^\s*WARNING:\s*\S", body):
        return True
    return False


def effective_risk_for_section(
    rev: Dict[str, Any],
    sec: Optional[Dict[str, Any]],
) -> str:
    """
    Stored ``risk`` when set (reviewer override). When missing/empty, infer from pipeline
    (high if WARNING / warning-kind, else low).
    """
    if rev:
        r = rev.get("risk")
        if r is not None and str(r).strip() != "":
            return _coerce_risk_level(r)
    if sec is not None:
        return "high" if _infer_pipeline_section_high_risk(sec) else "low"
    return "low"


def _risk_heading_colors(risk_level: str) -> Tuple[str, str]:
    """(border_hex, background_hex) for customer PDF / portal title bar."""
    if risk_level == "low":
        return "#15803d", "#f0fdf4"
    return "#be123c", "#fff1f2"


# Pipeline section titles (``detailed_sections[].title``) → customer-facing report titles.
# Keys: lowercase, whitespace-collapsed for robust matching.
_DARK_GENES_DISPLAY_TITLE: Dict[str, str] = {
    "hba analysis (alpha thalassemia - dosage)": "Alpha Thalassemia",
    "cyp21a2 analysis (cah - dosage)": "Congenital Adrenal Hyperplasia (CAH)",
    "expansion hunter (fragile x / fmr1)": "Fragile X",
    "dmd analysis (chrx:31.1m-33.3m)": "Duchenne Muscular Dystrophy (DMD)",
    "large svs (manta/gcnv - rest)": "Large Structural Variants and Copy Number Variants",
}


def dark_genes_display_title(raw: Optional[str]) -> str:
    """
    Map pipeline dark-gene section titles to report-facing names (PDF + portal).
    Unknown titles pass through unchanged.
    """
    t = (raw or "").strip()
    if not t:
        return "Section"
    if re.match(r"^smaca\s+check\b", t, re.I):
        return "Spinal Muscular Atrophy"
    key = " ".join(t.lower().split())
    return _DARK_GENES_DISPLAY_TITLE.get(key, t)


# Lines that are only decorative separators (e.g. runs of hyphens between blocks in *_detailed_report.txt)
_SECTION_SEPARATOR_LINE = re.compile(r"^[\s\-_–—]+$")


def _normalize_section_body(body: str) -> str:
    """Remove dash/underscore-only lines from section bodies (pipeline formatting noise)."""
    out: List[str] = []
    for line in (body or "").splitlines():
        if _SECTION_SEPARATOR_LINE.match(line):
            continue
        out.append(line)
    return "\n".join(out).strip()


# PDF: keep supplementary annex readable; raise if directors need more
_PDF_DETAILED_MAX_CHARS = 14_000

def _line_looks_like_summary_table_header(line: str) -> bool:
    """Header row from *_summary_report.txt (tab-, comma-, or pipe-separated)."""
    s = (line or "").strip()
    if not s:
        return False
    low = s.lower()
    if not (
        re.match(r"^sample\s*[\t,|]", low)
        or re.match(r"^sample\s{2,}\S", low)
    ):
        return False
    return any(
        x in low
        for x in ("paraphase", "smaca", "fragile", "hba", "cyp21", "large_sv", "qc_warn")
    )


def _line_looks_like_summary_data_row(line: str) -> bool:
    """Data row that belongs to the summary TSV (not SAMPLE: single-field lines)."""
    s = (line or "").strip()
    if not s:
        return False
    low = s.lower()
    if low.startswith("sample:") and "\t" not in s:
        return False
    if "\t" in s and len(s) > 12:
        return True
    if s.count(",") >= 4 and len(s) > 25:
        return True
    if s.count("|") >= 6 and len(s) > 50:
        return True
    return False


def _line_looks_like_summary_compact_pipeline_row(line: str) -> bool:
    """
    Single-line summary row (no separate header line) — some runs emit one dense row
    with many | fields and tool columns (Paraphase, SMAca, …).
    """
    s = (line or "").strip()
    if len(s) < 50:
        return False
    sl = s.lower()
    if not any(
        k in sl
        for k in (
            "paraphase",
            "smaca",
            "fragile",
            "smn1_cn",
            "smn2_cn",
            "fmr1",
            "hba",
            "cyp21",
            "large_sv",
            "qc_warn",
            "manta",
            "error:",
        )
    ):
        return False
    if s.count("|") >= 4:
        return True
    if "\t" in s and len(s) > 40:
        return True
    return False


def _strip_summary_tsv_from_text(text: str) -> str:
    """
    Remove summary TSV blocks (header + data row(s)) anywhere in the text.

    The same row as *_summary_report.txt is often pasted at the top of *_detailed_report.txt
    or embedded in the first «Overview» block; a prefix-only strip misses those cases.
    """
    if not (text or "").strip():
        return ""
    raw = (text or "").replace("\ufeff", "")
    lines = raw.splitlines()
    out: List[str] = []
    i = 0
    while i < len(lines):
        if _line_looks_like_summary_table_header(lines[i]):
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            while i < len(lines) and lines[i].strip() and _line_looks_like_summary_data_row(
                lines[i]
            ):
                i += 1
            continue
        if _line_looks_like_summary_compact_pipeline_row(lines[i]):
            i += 1
            continue
        out.append(lines[i])
        i += 1
    return "\n".join(out).strip()


def _detailed_text_for_dark_genes_pdf_parse(raw_detailed: str) -> str:
    """
    Text to feed ``parse_detailed_report_sections`` / synthetic sections.

    Prefer TSV-stripped text; if stripping removes everything (over-aggressive or file is
    mostly TSV-shaped lines), use a bounded copy of ``raw_detailed`` so the PDF path is
    not left with nothing to render.
    """
    raw = (raw_detailed or "").strip().replace("\ufeff", "")
    if not raw:
        return ""
    stripped = _strip_summary_tsv_from_text(raw)
    if stripped.strip():
        base = stripped
    else:
        base = raw
    if len(base) > _PDF_DETAILED_MAX_CHARS:
        return base[:_PDF_DETAILED_MAX_CHARS] + "\n\n[… truncated for PDF …]"
    return base


def _section_body_for_pdf(body: str) -> str:
    """Dash-only line cleanup + strip any pasted summary TSV from a section body."""
    return _normalize_section_body(_strip_summary_tsv_from_text(body or ""))


def sanitize_dark_genes_payload_for_pdf_render(report_data: Dict[str, Any]) -> None:
    """
    PDF render safety net: drop legacy ``report_summary`` / ``report_detailed`` from
    ``report.json`` when ``status != error`` so stale JSON cannot resurrect the raw TSV block.
    """
    dg = report_data.get("dark_genes")
    if not isinstance(dg, dict):
        return
    if dg.get("status") == "error":
        return
    dg.pop("report_summary", None)
    dg.pop("report_detailed", None)

def discover_dark_gene_visual_evidence(analysis_dir: str) -> Dict[str, Any]:
    """
    Find IGV ``igv-reports`` HTML and ExpansionHunter REViewer SVGs under the Nextflow outdir.

    Paths are **relative to analysis_dir** for use with ``/order/{id}/file/...`` (after adding
    ``analysis_dir`` to artifact roots). SMN1/SMN2/SMA loci live inside the unified HTML report.
    """
    base = (analysis_dir or "").strip()
    out: Dict[str, Any] = {
        "igv_report_html": None,
        "repeat_svgs": [],
        "snapshots_png": [],
    }
    if not base or not os.path.isdir(base):
        return out

    snap_roots = (
        os.path.join(base, "snapshots"),
        os.path.join(base, "results", "snapshots"),
    )
    igv_html: Optional[str] = None
    for snap in snap_roots:
        if not os.path.isdir(snap):
            continue
        htmls = sorted(glob.glob(os.path.join(snap, "*_visual_report.html")))
        if not htmls:
            htmls = sorted(glob.glob(os.path.join(snap, "*.html")))
        if htmls and not igv_html:
            igv_html = os.path.relpath(htmls[0], base)
        for png in sorted(glob.glob(os.path.join(snap, "*.png"))):
            rel = os.path.relpath(png, base)
            if rel not in out["snapshots_png"]:
                out["snapshots_png"].append(rel)
    out["igv_report_html"] = igv_html

    rep_roots = (
        os.path.join(base, "repeat"),
        os.path.join(base, "results", "repeat"),
    )
    for rep in rep_roots:
        if not os.path.isdir(rep):
            continue
        for svg in sorted(glob.glob(os.path.join(rep, "*.svg"))):
            rel = os.path.relpath(svg, base)
            if rel not in out["repeat_svgs"]:
                out["repeat_svgs"].append(rel)

    return out


def merge_visual_evidence_across_roots(search_roots: List[str]) -> Dict[str, Any]:
    """
    IGV HTML / repeat SVGs may live under a different directory than the first root that
    yielded ``summary/*.txt`` (e.g. Nextflow ``outdir`` vs mirrored ``output/``). Merge
    discoveries so relative paths like ``snapshots/foo.html`` resolve via
    ``/order/{id}/file/...`` (artifact roots include analysis + output).
    """
    out: Dict[str, Any] = {
        "igv_report_html": None,
        "repeat_svgs": [],
        "snapshots_png": [],
    }
    seen_svg: set = set()
    seen_png: set = set()
    for root in search_roots:
        if not root or not os.path.isdir(root):
            continue
        ev = discover_dark_gene_visual_evidence(root)
        if ev.get("igv_report_html") and not out["igv_report_html"]:
            out["igv_report_html"] = ev["igv_report_html"]
        for p in ev.get("repeat_svgs") or []:
            if p not in seen_svg:
                seen_svg.add(p)
                out["repeat_svgs"].append(p)
        for p in ev.get("snapshots_png") or []:
            if p not in seen_png:
                seen_png.add(p)
                out["snapshots_png"].append(p)
    return out


_MAX_READ_BYTES = 500_000
# Order matters: first pattern group with any match wins (prefer summary/ subfolder).
_SUMMARY_GLOBS = (
    "summary/*_summary_report.txt",
    "summary/*summary_report*.txt",
    "results/summary/*_summary_report.txt",
    "results/summary/*summary_report*.txt",
    "*_summary_report.txt",
    "*summary_report*.txt",
)
_DETAILED_GLOBS = (
    "summary/*_detailed_report.txt",
    "summary/*detailed_report*.txt",
    "results/summary/*_detailed_report.txt",
    "*_detailed_report.txt",
    "*detailed_report*.txt",
)
# SMAca process writes this (see docs/integrations/carrier_dark_genes.md); merged into detailed text.
_SMACA_CT_SIDECAR_GLOBS = (
    "summary/smaca_snp_counts.txt",
    "summary/SMAca/smaca_snp_counts.txt",
    "smaca/smaca_snp_counts.txt",
    "smaca_snp_counts.txt",
)


def _parse_smaca_ct_sidecar_text(raw: str) -> Optional[str]:
    """
    First usable line for portal/PDF: ``C_T=n,n`` or ``CT_counts=n,n``, or two integers.
    """
    if not (raw or "").strip():
        return None
    for line in raw.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if re.match(r"^C_T\s*=\s*\d+\s*,\s*\d+", s, re.I):
            return s
        if re.match(r"^CT_counts\s*=\s*\d+\s*,\s*\d+", s, re.I):
            return s
        if re.match(r"^SMN_C_T\s*=\s*\d+\s*,\s*\d+", s, re.I):
            return s
    lines0 = raw.strip().splitlines()
    if lines0:
        m = re.match(r"^(\d+)\s+(\d+)\s*$", lines0[0].strip())
        if m:
            return f"C_T={m.group(1)},{m.group(2)}"
    return None


def _inject_smaca_ct_line_into_detailed_text(detailed: str, ct_line: str) -> str:
    """
    Insert ``C_T=…`` into the SMAca block when the pipeline only wrote a sidecar file.
    Avoid duplicating if already present.
    """
    if not (detailed or "").strip() or not (ct_line or "").strip():
        return detailed
    if re.search(r"(?im)^\s*C_T\s*=", detailed):
        return detailed
    if re.search(r"(?im)^\s*CT_counts\s*=", detailed):
        return detailed
    inj = ct_line if ct_line.endswith("\n") else ct_line + "\n"
    # Prefer after C_Ratio= (same block as SMN1 cov fraction)
    m = re.search(r"(?im)^(\s*C_Ratio\s*=\s*[^\n]+\n)", detailed)
    if m:
        i = m.end()
        return detailed[:i] + inj + detailed[i:]
    # Else before SMN1_CN=
    m2 = re.search(r"(?im)^(\s*SMN1_CN\s*=\s*)", detailed)
    if m2:
        i = m2.start()
        return detailed[:i] + inj + detailed[i:]
    return detailed.rstrip() + "\n\n" + inj


def _read_smaca_ct_sidecar(base: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Returns (``C_T=…`` line or None, relative path of file read or None).
    """
    paths: List[str] = []
    for pat in _SMACA_CT_SIDECAR_GLOBS:
        paths.extend(glob.glob(os.path.join(base, pat)))
    if not paths:
        paths = _rglob_unique(base, ("**/smaca_snp_counts.txt",))
    if not paths:
        return None, None
    path = sorted(paths)[0]
    raw = _read_text_safe(path, 4096)
    line = _parse_smaca_ct_sidecar_text(raw)
    if not line:
        return None, None
    try:
        rel = os.path.relpath(path, base)
    except ValueError:
        rel = os.path.basename(path)
    return line, rel


def _read_text_safe(path: str, max_chars: int) -> str:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return f.read(max_chars + 1)[:max_chars]
    except OSError as e:
        logger.warning("[dark_genes] Could not read %s: %s", path, e)
        return ""


def _first_existing_glob(base: str, patterns: tuple) -> List[str]:
    """
    Try each pattern in order; return matches from the first pattern that yields files.
    (Avoid mixing summary/ and root files in one sort — prefer standard layout.)
    """
    for pat in patterns:
        batch: List[str] = []
        seen: set = set()
        for p in glob.glob(os.path.join(base, pat)):
            try:
                rp = os.path.realpath(p)
            except OSError:
                rp = p
            if rp in seen or not os.path.isfile(p):
                continue
            seen.add(rp)
            batch.append(p)
        if batch:
            return sorted(batch)
    return []


def _rglob_unique(base: str, patterns: Tuple[str, ...]) -> List[str]:
    """Recursive fallback when summary/ sits under an extra directory level."""
    out: List[str] = []
    seen = set()
    for pat in patterns:
        for p in glob.glob(os.path.join(base, pat), recursive=True):
            try:
                rp = os.path.realpath(p)
            except OSError:
                rp = p
            if rp in seen or not os.path.isfile(p):
                continue
            seen.add(rp)
            out.append(p)
    return sorted(out)


def parse_detailed_report_sections(text: str) -> List[Dict[str, Any]]:
    """
    Parse dark_gene_pipeline *_detailed_report.txt into sections for Review + PDF.

    Format (see modules/summary.nf): SAMPLE line, ===, !!! QUALITY WARNINGS !!!,
    then titled blocks ending with ':' (PARAPHASE RESULTS, SMAca CHECK, …).
    """
    raw = (text or "").strip()
    if not raw:
        return []

    lines = raw.splitlines()
    buf: List[str] = []
    title = "Overview"
    sections: List[Dict[str, Any]] = []

    def flush() -> None:
        nonlocal buf, title
        body = _normalize_section_body("\n".join(buf).strip())
        buf = []
        if not body:
            return
        tl = title.upper()
        kind = "normal"
        if "QUALITY" in tl and "WARNING" in tl:
            kind = "alert"
        elif "WARNING" in tl:
            kind = "warning"
        sections.append({"title": title.strip(), "body": body, "kind": kind})

    def is_banner(line: str) -> bool:
        s = line.strip()
        return bool(s) and s.startswith("!!!") and s.endswith("!!!") and len(s) < 160

    def is_section_header(line: str) -> bool:
        s = line.strip()
        if not s or len(s) > 130:
            return False
        if is_banner(line):
            return True
        if line.startswith("  ") or line.startswith("\t"):
            return False
        if not s.endswith(":"):
            return False
        # "KEY: value" on one line (e.g. SAMPLE: foo) — not a standalone header
        if re.match(r"^[A-Za-z0-9_]+\s*:\s*\S", s):
            return False
        return True

    for line in lines:
        if line.strip().startswith("=") and set(line.strip()) <= {"="}:
            continue
        if is_section_header(line):
            flush()
            st = line.strip()
            if is_banner(line):
                title = st.strip("! ").strip()
            else:
                title = st.rstrip(":").strip()
            continue
        buf.append(line)

    flush()
    return sections


def align_section_reviews(
    prev: Optional[List[Any]],
    n: int,
    sections: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """Align portal ``section_reviews`` to ``detailed_sections`` length (index-matched).

    ``risk`` defaults from the pipeline when missing: **low** unless the section has a
    ``WARNING:`` line or warning-kind title; reviewer-stored ``risk`` always wins.
    """
    out: List[Dict[str, Any]] = []
    for i in range(n):
        sec_i = sections[i] if sections and i < len(sections) else None
        if prev and i < len(prev) and isinstance(prev[i], dict):
            p = prev[i]
            pr = p.get("risk")
            if pr is not None and str(pr).strip() != "":
                risk = _coerce_risk_level(pr)
            else:
                risk = (
                    "high"
                    if (sec_i and _infer_pipeline_section_high_risk(sec_i))
                    else "low"
                )
            out.append(
                {
                    "approved": _coerce_approved_bool(p.get("approved")),
                    "notes": str(p.get("notes") or "")[:8000],
                    "risk": risk,
                }
            )
        else:
            risk = (
                "high" if (sec_i and _infer_pipeline_section_high_risk(sec_i)) else "low"
            )
            out.append({"approved": False, "notes": "", "risk": risk})
    return out


def merge_dark_genes_reviews(
    new_block: Dict[str, Any],
    prev_dg: Dict[str, Any],
) -> Dict[str, Any]:
    """Carry forward ``section_reviews`` after reprocess when section count still matches."""
    out = dict(new_block)
    prev_rev = prev_dg.get("section_reviews")
    if not isinstance(prev_rev, list):
        return out
    sections = out.get("detailed_sections")
    if not isinstance(sections, list) or len(sections) == 0:
        return out
    out["section_reviews"] = align_section_reviews(prev_rev, len(sections), sections)
    return out


def ensure_dark_genes_detailed_sections(result_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Fill ``dark_genes.detailed_sections`` from ``detailed_text`` when missing or empty.

    Older ``result.json`` / DB snapshots may only have ``detailed_text``; Review and PDF
    still need structured sections.
    Preserves ``section_reviews`` indices when re-parsing.
    """
    if not isinstance(result_data, dict):
        return result_data
    dg = result_data.get("dark_genes")
    if not isinstance(dg, dict):
        return result_data
    text = (dg.get("detailed_text") or "").strip()
    if not text:
        return result_data
    existing = dg.get("detailed_sections")
    if isinstance(existing, list) and len(existing) > 0:
        return result_data
    out = dict(result_data)
    dg2 = dict(dg)
    new_sections = parse_detailed_report_sections(text)
    dg2["detailed_sections"] = new_sections
    prev_rev = dg.get("section_reviews")
    if isinstance(prev_rev, list) and new_sections:
        dg2["section_reviews"] = align_section_reviews(prev_rev, len(new_sections), new_sections)
    out["dark_genes"] = dg2
    return out


def _infer_section_kind(sec: Dict[str, Any]) -> str:
    """Match parse_detailed_report_sections when ``kind`` is missing on stored sections."""
    k = sec.get("kind")
    if k in ("alert", "warning", "normal"):
        return str(k)
    tl = (sec.get("title") or "").strip().upper()
    if "QUALITY" in tl and "WARNING" in tl:
        return "alert"
    if "WARNING" in tl:
        return "warning"
    return "normal"


def _section_lab_review_only(sec: Dict[str, Any]) -> bool:
    """Pipeline QC / quality warning blocks — reviewers only; never customer PDF."""
    k = _infer_section_kind(sec)
    if k == "alert":
        return True
    t = (sec.get("title") or "").strip().lower()
    if "quality" in t and ("warning" in t or "warnings" in t):
        return True
    if "qc_warn" in t.replace(" ", ""):
        return True
    b = (sec.get("body") or "").strip().lower()
    # Same wording as summary QC_Warnings column / typical pipeline banner
    if len(b) < 800 and "median depth below" in b and "15x" in b:
        return True
    if len(b) < 800 and "alpha-cluster" in b and "sma" in b and "warning" in b:
        return True
    return False


def _section_is_entire_duplicate_summary_row(sec: Dict[str, Any]) -> bool:
    """One-line body that is the same dense summary row as *_summary_report.txt — omit from PDF."""
    body = (sec.get("body") or "").strip()
    if not body or "\n" in body or "\r" in body:
        return False
    return _line_looks_like_summary_compact_pipeline_row(body)


def _section_review_at(
    section_reviews: Optional[List[Any]], i: int
) -> Dict[str, Any]:
    """Index ``i`` review dict, or ``{}`` if missing / not a dict (matches gate + render)."""
    if not section_reviews or i < 0 or i >= len(section_reviews):
        return {}
    x = section_reviews[i]
    return x if isinstance(x, dict) else {}


def _dl_kv_pdf_rows(rows: List[Tuple[str, str]]) -> str:
    """Portal-style label/value rows for WeasyPrint (inline styles; no portal CSS)."""
    parts: List[str] = []
    for lab, val in rows:
        parts.append(
            "<tr>"
            f'<td style="padding:3px 12px 3px 0;font-size:8pt;color:#64748b;'
            f'vertical-align:top;white-space:nowrap;">{html.escape(lab)}</td>'
            f'<td style="padding:3px 0;font-size:8pt;color:#0f172a;">{html.escape(val)}</td>'
            "</tr>"
        )
    return (
        '<table style="margin:4px 0 0;border-collapse:collapse;width:100%;">'
        f'{"".join(parts)}</table>'
    )


def _smaca_extract_snp_ct_counts(raw: str) -> Optional[Tuple[str, str]]:
    """
    Best-effort SMN1/SMN2-discriminating SNP read counts from pipeline text.

    Tries several ``key=value`` styles (and a loose g.27134 depth line) so the
    portal/PDF can show C/T next to ``C_Ratio`` when the detailed report includes them.
    """
    if not (raw or "").strip():
        return None
    patterns = (
        r"CT_counts\s*=\s*(\d+)\s*,\s*(\d+)",
        r"C_T\s*=\s*(\d+)\s*,\s*(\d+)",
        r"SMN_C_T\s*=\s*(\d+)\s*,\s*(\d+)",
    )
    for pat in patterns:
        m = re.search(pat, raw, re.I)
        if m:
            return (m.group(1), m.group(2))
    mc = re.search(r"(?im)^\s*C_reads\s*=\s*(\d+)\s*$", raw)
    mt = re.search(r"(?im)^\s*T_reads\s*=\s*(\d+)\s*$", raw)
    if mc and mt:
        return (mc.group(1), mt.group(1))
    m = re.search(
        r"27134[^\n]{0,220}?(?:AD|DP|depth|counts?)\s*[=:]\s*(\d+)\s*[,;/]\s*(\d+)",
        raw,
        re.I,
    )
    if m:
        return (m.group(1), m.group(2))
    return None


def _try_smaca_kv_html(title: str, body: str) -> Optional[str]:
    """Match portal ``tryRenderSmacaCheckSection`` (SMN1/SMN2/Silent Carrier/Ratio)."""
    if not re.search(r"SMAca\s*CHECK", title, re.I):
        return None
    raw = body or ""
    m1e = re.search(r"SMN1_CN_est\s*=\s*(\S+)", raw, re.I)
    m2e = re.search(r"SMN2_CN_est\s*=\s*(\S+)", raw, re.I)
    m1 = re.search(r"SMN1_CN\s*=\s*(\S+)", raw, re.I)
    m2 = re.search(r"SMN2_CN\s*=\s*(\S+)", raw, re.I)
    m_silent = re.search(r"SilentCarrier\s*=\s*([^\s(]+)(?:\s*\([^)]*\))?", raw, re.I)
    # Pipeline: C_Ratio = avg_cov_SMN1 / (avg_cov_SMN1 + avg_cov_SMN2), not allele C/T.
    m_cr = re.search(r"C_Ratio\s*=\s*(\S+)", raw, re.I)
    m_cov = re.search(r"Cov\s*\(\s*1\s*,\s*2\s*\)\s*=\s*(.+)", raw, re.I)
    # Optional: g.27134-style SNP read split (emit from pipeline as key=value).
    m_ct_ratio = re.search(r"CT_Ratio\s*=\s*(\S+)", raw, re.I)
    ct_pair = _smaca_extract_snp_ct_counts(raw)
    if not m1e and not m2e and not m1 and not m2 and not m_silent and not m_cr and not m_cov and not m_ct_ratio and not ct_pair:
        return None
    rows: List[Tuple[str, str]] = []
    if m1e:
        rows.append(("SMN1 CNV (est.)", m1e.group(1)))
    elif m1:
        rows.append(("SMN1 CNV", m1.group(1)))
    if m2e:
        rows.append(("SMN2 CNV (est.)", m2e.group(1)))
    elif m2:
        rows.append(("SMN2 CNV", m2.group(1)))
    if m_silent:
        rows.append(("Silent Carrier", m_silent.group(1)))
    if m_cr:
        rows.append(("SMN1 cov fraction", m_cr.group(1)))
    if ct_pair:
        rows.append(("SNP C,T counts", f"{ct_pair[0]} / {ct_pair[1]}"))
    if m_ct_ratio:
        rows.append(("SNP C/T ratio", m_ct_ratio.group(1)))
    if m_cov:
        rows.append(("Cov(1,2)", re.sub(r"\s+", " ", m_cov.group(1).strip())))
    return _dl_kv_pdf_rows(rows)


def _dosage_title_matches(title: str) -> bool:
    """Match portal ``dosageAnalysisTitleMatches``."""
    u = (title or "").strip()
    if not u:
        return False
    if re.search(r"HBA\s+ANALYSIS", u, re.I):
        return True
    if re.search(r"CYP21A2\s+ANALYSIS", u, re.I):
        return True
    if re.search(r"Alpha\s*Thalassemia", u, re.I) and re.search(r"Dosage", u, re.I):
        return True
    if re.search(r"CYP21A2", u, re.I) and re.search(r"CAH", u, re.I):
        return True
    if re.search(r"CAH", u, re.I) and re.search(r"Dosage", u, re.I) and re.search(r"CYP21", u, re.I):
        return True
    return False


def _try_dosage_kv_html(title: str, body: str) -> Optional[str]:
    """Match portal ``tryRenderDosageAnalysisSection`` (Estimated CNV, Ratio, Warning)."""
    if not _dosage_title_matches(title):
        return None
    raw = body or ""
    m_est = re.search(r"^\s*Est_CN\s*=\s*(\S+)", raw, re.I | re.M)
    m_ratio = re.search(r"^\s*Ratio\s*=\s*(\S+)", raw, re.I | re.M)
    m_warn = re.search(r"^\s*WARNING:\s*(.+)", raw, re.I | re.M)
    if not m_est and not m_ratio and not m_warn:
        return None
    rows: List[Tuple[str, str]] = []
    if m_est:
        rows.append(("Estimated CNV", m_est.group(1)))
    if m_ratio:
        rows.append(("Ratio", m_ratio.group(1)))
    if m_warn:
        rows.append(("Warning", re.sub(r"\s+", " ", m_warn.group(1).strip())))
    return _dl_kv_pdf_rows(rows)


def _try_generic_pipeline_kv_html(body: str) -> Optional[str]:
    """
    Fallback: KEY=value lines, WARNING:, or single-line prose / ``Gene:value`` (e.g. FMR1:13/13).
    Keeps PDF readable without raw monospaced pipeline dumps when SMA/dosage matchers miss.
    """
    raw = (body or "").strip()
    if not raw:
        return None
    lines = [
        ln.strip()
        for ln in raw.splitlines()
        if ln.strip() and not re.match(r"^[\s\-_–—]+$", ln)
    ]
    if not lines:
        return None
    rows: List[Tuple[str, str]] = []
    leftover: List[str] = []
    for ln in lines:
        if re.match(r"^WARNING:\s*", ln, re.I):
            rows.append(("Warning", re.sub(r"^WARNING:\s*", "", ln, flags=re.I).strip()))
            continue
        m_eq = re.match(r"^([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+)$", ln)
        if m_eq:
            key, val = m_eq.group(1), m_eq.group(2).strip()
            rows.append((key.replace("_", " "), val))
            continue
        m_colon = re.match(r"^([A-Za-z0-9]+):(.+)$", ln)
        if m_colon and len(m_colon.group(1)) <= 12:
            rows.append((m_colon.group(1), m_colon.group(2).strip()))
            continue
        leftover.append(ln)
    if rows and not leftover:
        return _dl_kv_pdf_rows(rows)
    if rows and leftover:
        return _dl_kv_pdf_rows(rows) + (
            '<p style="white-space:pre-wrap;font-size:8pt;margin:6px 0 0;line-height:1.45;'
            f'color:#334155;">{html.escape(chr(10).join(leftover))}</p>'
        )
    if len(lines) == 1:
        return (
            f'<p style="font-size:8pt;color:#334155;margin:8px 0 0;line-height:1.45;">'
            f"{html.escape(lines[0])}</p>"
        )
    return (
        '<p style="white-space:pre-wrap;font-size:8pt;margin:8px 0 0;line-height:1.45;'
        f'color:#334155;">{html.escape(chr(10).join(lines))}</p>'
    )


def _section_body_portal_html_for_pdf(sec: Dict[str, Any]) -> str:
    """
    Same lab-facing labels as the portal Review tab (``tryRenderSmacaCheckSection`` /
    ``tryRenderDosageAnalysisSection``), not raw ``*_detailed_report.txt`` monospaced lines.
    """
    title = (sec.get("title") or "").strip()
    body = _section_body_for_pdf(sec.get("body") or "")
    h = _try_smaca_kv_html(title, body)
    if h:
        return h
    h = _try_dosage_kv_html(title, body)
    if h:
        return h
    gh = _try_generic_pipeline_kv_html(body)
    return gh or ""


def _dark_genes_finding_card_html(
    title_escaped: str,
    analysis_inner_html: str,
    notes_plain: str,
    *,
    margin_top: str = "25px",
    risk_level: str = "high",
) -> str:
    """
    Match ``carrier_EN.html`` Detailed Interpretations: left accent bar (green if ``risk_level``
    is ``low``, red if ``high`` — reviewer-set), uppercase sub-labels, dashed separator.
    """
    border, bg = _risk_heading_colors(_coerce_risk_level(risk_level))
    head = (
        f'<div style="font-size:11pt;font-weight:500;color:#0f172a;margin-bottom:15px;'
        f'padding:4px 0 4px 12px;border-left:6px solid {border};line-height:1.2;'
        f'background-color:{bg};display:block;">{title_escaped}</div>'
    )
    analysis_wrap = ""
    inner = (analysis_inner_html or "").strip()
    if inner:
        sub_test = (
            '<div style="font-size:8pt;font-weight:500;color:#64748b;text-transform:uppercase;'
            'letter-spacing:0.02em;margin-top:10px;">Test analysis</div>'
        )
        analysis_wrap = (
            f'<div style="margin-bottom:12px;">{sub_test}'
            f'<div style="margin-top:2px;">{analysis_inner_html}</div></div>'
        )
    notes_block = ""
    notes = (notes_plain or "").strip()
    if notes:
        sub_int = (
            '<div style="font-size:8pt;font-weight:500;color:#64748b;text-transform:uppercase;'
            'letter-spacing:0.02em;margin-top:10px;">Interpretation</div>'
        )
        esc = html.escape(notes, quote=True)
        notes_block = (
            f'<div>{sub_int}'
            f'<p style="font-size:8pt;color:#334155;margin:2px 0 0 0;text-align:justify;'
            f'line-height:1.3;white-space:pre-wrap;">{esc}</p></div>'
        )
    return (
        f'<div style="margin-top:{margin_top};padding:10px 0;page-break-inside:avoid;'
        'border-bottom:1px dashed #e2e8f0;">'
        f"{head}{analysis_wrap}{notes_block}</div>"
    )


def detailed_sections_to_pdf_html(
    sections: List[Dict[str, Any]],
    section_reviews: Optional[List[Dict[str, Any]]] = None,
    *,
    filter_by_approval: bool = False,
) -> str:
    """
    Escaped HTML blocks for WeasyPrint (carrier_*.html uses |safe).

    Quality warning blocks (``kind: alert``) are omitted — lab review only.
    The **Overview** block is omitted from the customer PDF (keep it in the portal for audit only).
    One-line duplicate summary rows are omitted.

    When ``filter_by_approval=True`` (customer PDF), only sections with
    ``section_reviews[i].approved`` are included. Title accent is **green** if
    effective risk is ``low``, **red** if ``high`` (reviewer override, or pipeline
    ``WARNING:`` / warning-kind section). Each block
    matches **Detailed Interpretations** layout (**Test analysis** + **Interpretation**
    sub-labels, dashed separators).
    When ``filter_by_approval=False``, every eligible section includes full body text and
    optional **Notes** lines beneath.
    """
    if not sections:
        return ""
    parts: List[str] = []
    first_approved_card = True
    for i, sec in enumerate(sections):
        if _section_lab_review_only(sec):
            continue
        if _section_is_entire_duplicate_summary_row(sec):
            continue
        title_l = (sec.get("title") or "").strip().lower()
        if title_l == "overview":
            continue
        rev = _section_review_at(section_reviews, i)
        if filter_by_approval and not _coerce_approved_bool(rev.get("approved")):
            continue
        notes = (rev.get("notes") or "").strip()

        t = html.escape(dark_genes_display_title(sec.get("title")), quote=True)
        body_for_pdf = _section_body_for_pdf(sec.get("body") or "")
        b = html.escape(body_for_pdf, quote=True)
        kind = _infer_section_kind(sec)
        if kind == "warning":
            box = (
                "margin-top:10px;border-radius:6px;border:1px solid #f59e0b;"
                "background:#fffbeb;padding:10px 12px;"
            )
            tit = "font-weight:700;font-size:9pt;color:#b45309;margin:0 0 6px;"
        else:
            box = (
                "margin-top:10px;border-radius:6px;border:1px solid #cbd5e1;"
                "background:#f8fafc;padding:10px 12px;"
            )
            tit = "font-weight:700;font-size:9pt;color:#0f172a;margin:0 0 6px;"

        if filter_by_approval:
            # Customer PDF: same layout as Detailed Interpretations (condition-heading + sub-labels).
            portal_block = _section_body_portal_html_for_pdf(sec)
            mt = "12px" if first_approved_card else "25px"
            first_approved_card = False
            parts.append(
                _dark_genes_finding_card_html(
                    t,
                    portal_block,
                    notes,
                    margin_top=mt,
                    risk_level=effective_risk_for_section(rev, sec),
                )
            )
            continue

        notes_html = ""
        if notes:
            notes_html = (
                f'<p style="font-size:8pt;color:#475569;margin:8px 0 0;line-height:1.4;">'
                f"<strong>Notes:</strong> {html.escape(notes, quote=True)}</p>"
            )
        parts.append(
            f'<div style="{box}">'
            f'<div style="{tit}">{t}</div>'
            f'<pre style="white-space:pre-wrap;font-size:8pt;line-height:1.45;'
            f"margin:0;color:#334155;font-family:ui-monospace,Consolas,monospace\">{b}</pre>"
            f"{notes_html}</div>"
        )
    return "\n".join(parts)


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
    if not summary_paths:
        summary_paths = _rglob_unique(
            base,
            ("**/*_summary_report.txt", "**/*summary_report*.txt"),
        )
    detailed_paths = _first_existing_glob(base, _DETAILED_GLOBS)
    if not detailed_paths:
        detailed_paths = _rglob_unique(
            base,
            ("**/*_detailed_report.txt", "**/*detailed_report*.txt"),
        )

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

    smaca_ct_line, smaca_ct_sidecar = _read_smaca_ct_sidecar(base)
    if smaca_ct_line and detailed_excerpt:
        detailed_excerpt = _inject_smaca_ct_line_into_detailed_text(
            detailed_excerpt, smaca_ct_line
        )

    status = "found" if summary_text else "partial"

    detailed_sections: List[Dict[str, Any]] = []
    if detailed_excerpt:
        try:
            detailed_sections = parse_detailed_report_sections(detailed_excerpt)
        except Exception as e:
            logger.warning("[dark_genes] parse_detailed_report_sections: %s", e)

    return {
        "status": status,
        "sample_name": (sample_name or "").strip(),
        "summary_file": os.path.basename(summary_paths[0]) if summary_paths else None,
        "detailed_file": os.path.basename(detailed_paths[0]) if detailed_paths else None,
        "summary_paths": [os.path.relpath(p, base) for p in summary_paths],
        "detailed_paths": [os.path.relpath(p, base) for p in detailed_paths],
        "summary_text": summary_text,
        "detailed_text": detailed_excerpt,
        "detailed_sections": detailed_sections,
        "visual_evidence": discover_dark_gene_visual_evidence(base),
        **(
            {"smaca_ct_sidecar": smaca_ct_sidecar}
            if smaca_ct_sidecar
            else {}
        ),
    }


def dark_genes_for_pdf(report_block: Dict[str, Any]) -> Dict[str, Any]:
    """
    Subset for customer PDF (written into ``report.json`` — not the full ``dark_genes`` blob).

    Reads ``detailed_text``, ``detailed_sections``, and ``section_reviews`` (index-aligned).
    The customer PDF **only includes sections marked approved** in the portal; unapproved
    sections are omitted (Overview / QC-only / duplicate summary rows are never included).
    For each approved section, ``report_detailed_html`` contains the **portal section title**,
    **portal-style parsed content** from the pipeline body (KV / labels like the Review tab),
    then **reviewer notes** — not a raw monospaced pipeline block.

    Only ``report_detailed_html`` is emitted. ``status == error`` returns a short ``report_summary``.
    """
    if not report_block or report_block.get("status") == "not_found":
        return {}
    if report_block.get("status") == "error":
        msg = (report_block.get("message") or "unknown error")[:4000]
        return {
            "status": "error",
            "report_summary": f"Supplementary dark-gene analysis metadata could not be loaded: {msg}",
        }
    raw_detailed = (report_block.get("detailed_text") or "").strip().replace("\ufeff", "")
    detailed = _detailed_text_for_dark_genes_pdf_parse(raw_detailed)

    parsed = parse_detailed_report_sections(detailed) if detailed else []
    section_reviews: Optional[List[Dict[str, Any]]] = None
    if isinstance(report_block, dict) and "section_reviews" in report_block:
        raw = report_block.get("section_reviews")
        section_reviews = []
        for x in (raw if isinstance(raw, list) else []):
            if isinstance(x, dict):
                item: Dict[str, Any] = {
                    "approved": _coerce_approved_bool(x.get("approved")),
                    "notes": str(x.get("notes") or "")[:8000],
                }
                xr = x.get("risk")
                if xr is not None and str(xr).strip() != "":
                    item["risk"] = _coerce_risk_level(xr)
                section_reviews.append(item)
            else:
                section_reviews.append({"approved": False, "notes": ""})

    stored = report_block.get("detailed_sections")
    # Always prefer on-disk ``detailed_sections`` (same indices as portal ``section_reviews``).
    # Falling back to re-parse from ``detailed_text`` can change section count/order and mis-apply approvals.
    if isinstance(stored, list) and len(stored) > 0:
        sections = stored
    else:
        sections = parsed

    # Parser expects titled blocks (lines ending with ":"). Some runs are free-form text only.
    if not sections and detailed:
        sections = [
            {
                "title": "Supplementary detail",
                "body": detailed,
                "kind": "normal",
            }
        ]

    if (
        section_reviews is None
        or not isinstance(section_reviews, list)
        or len(section_reviews) != len(sections)
    ):
        section_reviews = align_section_reviews(
            section_reviews if isinstance(section_reviews, list) else None,
            len(sections),
            sections,
        )

    report_detailed_html = (
        detailed_sections_to_pdf_html(sections, section_reviews, filter_by_approval=True)
        if sections
        else ""
    )

    if not (report_detailed_html or "").strip():
        if sections and isinstance(section_reviews, list) and len(section_reviews) == len(sections):
            logger.info(
                "[dark_genes_for_pdf] no approved sections for customer PDF (detailed_sections=%d)",
                len(sections),
            )
        else:
            logger.warning(
                "[dark_genes_for_pdf] empty supplemental HTML (sections=%d, "
                "section_reviews_len=%s, status=%s)",
                len(sections),
                len(section_reviews) if isinstance(section_reviews, list) else None,
                report_block.get("status"),
            )
        return {}
    return {
        "status": report_block.get("status"),
        "report_detailed_html": report_detailed_html,
    }
