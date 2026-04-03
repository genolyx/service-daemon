"""
Report Generator Module

리뷰어가 확정한 결과를 바탕으로 최종 리포트를 생성합니다.
Carrier_result/generate.py의 Jinja2 + WeasyPrint 패턴을 통합하여
다국어 PDF 리포트를 생성합니다.

생성 파일:
    - report.json       : 최종 확정 결과 (Portal → 고객 전달용)
    - report.html       : HTML 리포트 (미리보기용)
    - report.pdf        : PDF 리포트 (다국어 지원)

report.json 구조:
    - report_metadata: order_id, report_date, hospital, doctor, language, …
    - primary_patient: 환자 정보 + findings[] (Jinja carrier_*.html — 템플릿이 변이 표시에 사용)
    - partner: 파트너 정보 + findings[] (couples 템플릿)
    - findings: couples 전용, primary+partner 병합 (carrier_couples_*.html 상세 해석)
    - genes_evaluated_count
    - carrier_status, confirmed_variants, disease_groups, qc_summary, reviewer
    - dark_genes (optional): report_detailed_html (+ error-only report_summary) from result.json for PDF (approved detailed sections only; Overview/QC blocks omitted)
"""

import html
import os
import re
import json
import logging
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple
from zoneinfo import ZoneInfo

from ...datetime_kst import now_kst_date_iso, now_kst_iso

logger = logging.getLogger(__name__)


def _dedupe_result_json_paths(paths: Sequence[Optional[str]]) -> List[str]:
    """Absolute, deduplicated candidate paths for result.json (order preserved)."""
    seen: set = set()
    out: List[str] = []
    for p in paths:
        if not p:
            continue
        ap = os.path.normpath(os.path.abspath(p))
        if ap in seen:
            continue
        seen.add(ap)
        out.append(ap)
    return out


def _dark_genes_result_json_paths_for_pdf(
    primary_result_json: str,
    extra_result_json_paths: Optional[Sequence[str]] = None,
) -> List[str]:
    """
    ``result.json`` next to ``report.json`` (primary) **first** — portal saves and
    ``section_reviews`` usually live there. Remaining candidates follow, **newest mtime first**,
    so we still pick up fresh pipeline text when the primary path lags a duplicate copy.
    """
    ordered = _dedupe_result_json_paths(
        [primary_result_json] + list(extra_result_json_paths or [])
    )
    existing = [p for p in ordered if os.path.isfile(p)]
    primary_abs = os.path.normpath(os.path.abspath(primary_result_json))
    primary_file = None
    for p in existing:
        if os.path.normpath(os.path.abspath(p)) == primary_abs:
            primary_file = p
            break
    others = [p for p in existing if p != primary_file]
    try:
        others.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    except OSError:
        pass
    if primary_file:
        return [primary_file] + others
    return others


def _build_dark_genes_for_pdf_from_result_json_candidates(
    ordered_paths: List[str],
    log_context: str,
) -> Tuple[Optional[Dict[str, Any]], Optional[str], bool]:
    """
    Build customer-PDF ``dark_genes`` via ``dark_genes_for_pdf`` from one or more
    ``result.json`` paths.

    **Why not "first path that returns non-empty"?** Report-adjacent ``result.json`` may be
    older than a pipeline copy but sort first; or the newest file may have fresh
    ``detailed_sections`` while portal ``section_reviews`` only exist on another path.
    We take **detailed_text / detailed_sections** from the **newest mtime** usable blob, and
    **section_reviews** from the **newest mtime** file that actually contains that key
    (portal saves), then run ``dark_genes_for_pdf`` once on the merged dict.

    Returns ``(pdf_dg, source_description, any_file_opened)``.
    """
    from .dark_genes import dark_genes_for_pdf

    rows: List[Tuple[float, str, Dict[str, Any]]] = []
    any_file_opened = False
    for rp in ordered_paths:
        if not os.path.isfile(rp):
            continue
        any_file_opened = True
        try:
            with open(rp, "r", encoding="utf-8") as f:
                rd = json.load(f)
        except Exception as e:
            logger.warning("[%s] could not read %s: %s", log_context, rp, e)
            continue
        if not isinstance(rd, dict):
            continue
        dg = rd.get("dark_genes")
        if not isinstance(dg, dict) or dg.get("status") == "not_found":
            continue
        try:
            mt = os.path.getmtime(rp)
        except OSError:
            mt = 0.0
        rows.append((mt, rp, dg))

    if not rows:
        return None, None, any_file_opened

    rows.sort(key=lambda x: x[0], reverse=True)
    _, base_path, base_dg = rows[0]

    # Portal review state: prefer newest file among those that include section_reviews.
    with_reviews = [(mt, p, dg) for mt, p, dg in rows if "section_reviews" in dg]
    if with_reviews:
        with_reviews.sort(key=lambda x: x[0], reverse=True)
        # Prefer non-empty review lists when multiple files have the key (stale [] vs real saves).
        non_empty = [r for r in with_reviews if isinstance(r[2].get("section_reviews"), list) and len(r[2]["section_reviews"]) > 0]
        _mt, rev_path, rev_dg = (non_empty[0] if non_empty else with_reviews[0])
        merged = dict(base_dg)
        merged["section_reviews"] = rev_dg.get("section_reviews")
        pdf_dg = dark_genes_for_pdf(merged)
        if pdf_dg:
            src = (
                f"{base_path} + section_reviews:{rev_path}"
                if os.path.normpath(base_path) != os.path.normpath(rev_path)
                else base_path
            )
            return pdf_dg, src, any_file_opened
        logger.warning(
            "[%s] merged dark_genes (base=%s reviews=%s) produced empty PDF HTML",
            log_context,
            base_path,
            rev_path,
        )

    pdf_dg = dark_genes_for_pdf(base_dg)
    if pdf_dg:
        return pdf_dg, base_path, any_file_opened

    logger.warning(
        "[%s] dark_genes_for_pdf returned empty after trying merge + base (%s)",
        log_context,
        base_path,
    )
    return None, None, any_file_opened


def _merge_dark_genes_from_result_json_for_pdf(
    report_data: Dict[str, Any],
    report_json_path: str,
    extra_result_json_paths: Optional[Sequence[str]] = None,
) -> None:
    """
    Refresh ``report_data[\"dark_genes\"]`` from ``result.json`` next to ``report.json``.

    PDF rendering only reads ``report.json``. If that snapshot omitted the supplement (stale
    generate, old daemon, or hand-edited file), WeasyPrint would still omit the section — so we
    always re-merge from ``result.json`` at PDF time when present.

    When ``CARRIER_SCREENING_REPORT_OUTPUT_ROOT`` is set, ``report.json`` may live under that
    tree while the pipeline still wrote ``result.json`` only under ``job.output_dir``. Pass
    ``extra_result_json_paths`` (e.g. the work-tree ``result.json``) so the supplement still loads.
    """
    primary = os.path.join(os.path.dirname(os.path.abspath(report_json_path)), "result.json")
    deduped_all = _dedupe_result_json_paths(
        [primary] + list(extra_result_json_paths or [])
    )
    candidates = _dark_genes_result_json_paths_for_pdf(
        primary, extra_result_json_paths
    )
    pdf_dg, src, any_file = _build_dark_genes_for_pdf_from_result_json_candidates(
        candidates, "generate_report_pdf"
    )
    if pdf_dg and src:
        report_data["dark_genes"] = pdf_dg
        hlen = len((pdf_dg.get("report_detailed_html") or ""))
        logger.info(
            "[generate_report_pdf] merged dark_genes from %s (report_detailed_html chars=%s)",
            src,
            hlen,
        )
        return

    report_data.pop("dark_genes", None)
    if any_file:
        logger.warning(
            "[generate_report_pdf] no usable dark_genes supplement after trying (newest-first): %s",
            candidates,
        )
    elif deduped_all:
        logger.warning(
            "[generate_report_pdf] no result.json on disk for dark_genes merge (looked for %s)",
            deduped_all,
        )


# ══════════════════════════════════════════════════════════════
# Template selection (order params → carrier_{lang}.html vs carrier_couples_{lang}.html)
# ══════════════════════════════════════════════════════════════

def _carrier_order_flat(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    포털은 임상 메타를 params.carrier 에 저장한다.
    템플릿 판별은 carrier 키를 최상위와 병합해 읽는다.
    """
    if not params:
        return {}
    c = params.get("carrier")
    if isinstance(c, dict) and c:
        merged = dict(params)
        merged.update(c)
        return merged
    return dict(params)


def carrier_report_template_kind(params: Dict[str, Any]) -> Optional[str]:
    """
    PDF 템플릿 종류.
    - 'standard' → carrier_<lang>.html (Primary: Carrier screening standard)
    - 'couples'  → carrier_couples_<lang>.html (Other + CouplesCarrier)
    - None       → 위 둘 외 (Exome 등) — PDF 템플릿 미지원
    """
    if not params:
        return None
    params = _carrier_order_flat(params)
    tc = str(params.get("test_category") or "").strip()
    ot = str(params.get("other_test_type") or "").strip()
    pc = str(params.get("package_code") or "").strip()

    # Whole Exome orders (portal: service_code whole_exome, package_code WholeExome)
    if pc == "WholeExome":
        return "standard"

    if tc == "standard_carrier":
        return "standard"
    if tc == "other" and ot == "CouplesCarrier":
        return "couples"
    if pc == "CouplesCarrier":
        return "couples"
    if pc == "CarrierScreening" and tc != "other":
        return "standard"
    if tc == "other":
        return None
    if pc == "CarrierScreening":
        return "standard"
    return None


def report_languages_from_order(params: Dict[str, Any]) -> Optional[List[str]]:
    """
    order.params.report_language → 템플릿에 있는 EN/CN/KO 만 허용.
    ID / Other 등은 템플릿이 없으므로 None.
    """
    if not params:
        return None
    params = _carrier_order_flat(params)
    rl = str(params.get("report_language") or "").strip().upper()
    if rl in ("EN", "CN", "KO"):
        return [rl]
    return None


# ══════════════════════════════════════════════════════════════
# Jinja template JSON (Carrier_result / data/carrier_report/*.html)
# Templates iterate data.primary_patient.findings (not confirmed_variants).
# ══════════════════════════════════════════════════════════════

def _report_date_human_kst() -> str:
    """Matches sample JSON like 'Dec 24, 2025' for report cover."""
    return datetime.now(ZoneInfo("Asia/Seoul")).strftime("%b %d, %Y")


# ClinVar CLNDN often lists multiple names joined by |; also contains not_specified / not_provided.
_CLINVAR_DISORDER_SKIP = frozenset(
    {
        ".",
        "not_provided",
        "not_specified",
        "not_available",
        "see_cases",
        "inborn_genetic_diseases",
        "inborn_genetic_disease",
        "no_classification_for_the_single_variant",
        "notspecified",
        "notprovided",
    }
)


def _normalize_disorder_display(name: str) -> str:
    """Underscores → spaces for titles like Retinitis_pigmentosa_39."""
    s = (name or "").strip().replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def _best_disorder_from_clinvar_field(raw: str) -> str:
    """
    Pick a single human-readable condition from ClinVar/HGMD disease text.
    Splits on | and ;, drops placeholder tokens, returns the first substantive label.
    """
    if not raw or not str(raw).strip():
        return ""
    parts = re.split(r"[|;]", str(raw))
    candidates: List[str] = []
    for p in parts:
        t = p.strip()
        if not t:
            continue
        tl = t.lower().replace("_", "")
        if tl in _CLINVAR_DISORDER_SKIP:
            continue
        if tl.startswith("not_") and ("specified" in tl or "provided" in tl):
            continue
        if len(tl) < 2:
            continue
        candidates.append(t)

    if not candidates:
        return ""

    # Longest token often matches the specific syndrome vs a shorter alias in the same field
    pick = max(candidates, key=len)
    return _normalize_disorder_display(pick)


def _finalize_disorder_label(label: str) -> str:
    """
    One line for the CONDITION column. If the source stuffed a full ClinVar CLNDN
    (pipe-joined) into a disease name, split and pick a single label.
    """
    label = (label or "").strip()
    if not label:
        return ""
    if "|" in label or ";" in label:
        one = _best_disorder_from_clinvar_field(label)
        if one:
            return one
    return _normalize_disorder_display(label)


def _coalesce_inheritance_for_report_variant(v: Dict[str, Any]) -> str:
    """Same logic as review result.json: top-level or first diseases[].inheritance."""
    inh = (v.get("inheritance") or "").strip()
    if inh:
        return inh
    for d in v.get("diseases") or []:
        if isinstance(d, dict):
            x = (d.get("inheritance") or "").strip()
            if x:
                return x
    return ""


def _normalize_inheritance_lookup_key(label: str) -> str:
    s = (label or "").strip().lower().replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def _load_panel_disease_inheritance_by_name(panel_json_path: str) -> Dict[str, str]:
    """
    disease_gene_mapping.json diseases[].disease_name → inheritance (AR/AD/…).
    Used when variant has a clean disorder title but no inheritance on the object.
    """
    out: Dict[str, str] = {}
    try:
        with open(panel_json_path, "r", encoding="utf-8") as f:
            raw = json.load(f)
    except OSError:
        return out
    for row in raw.get("diseases") or []:
        if not isinstance(row, dict):
            continue
        inh = (row.get("inheritance") or "").strip()
        if not inh:
            continue
        for key in (
            row.get("disease_name"),
            row.get("name"),
            row.get("disease_name_ko"),
            row.get("disease_name_cn"),
        ):
            if not key or not str(key).strip():
                continue
            k = _normalize_inheritance_lookup_key(str(key))
            if k:
                out[k] = inh
    return out


def _inheritance_from_gene_panel_mapper(gene: str, mapper: Any) -> str:
    if not mapper or not (gene or "").strip():
        return ""
    info = mapper.lookup((gene or "").strip().upper())
    diseases = info.get("diseases") or []
    if diseases and isinstance(diseases[0], dict):
        return (diseases[0].get("inheritance") or "").strip()
    return ""


def _inheritance_from_panel_disease_name_map(
    disorder_label: str, name_map: Dict[str, str]
) -> str:
    """Exact match, then longest panel disease_name substring contained in disorder_label."""
    if not disorder_label or disorder_label == "Unknown disorder" or not name_map:
        return ""
    k = _normalize_inheritance_lookup_key(disorder_label)
    if k in name_map:
        return name_map[k]
    best_len = 0
    best_inh = ""
    for dk, inh in name_map.items():
        if len(dk) < 8:
            continue
        if dk in k and len(dk) > best_len:
            best_len = len(dk)
            best_inh = inh
    return best_inh


def _resolve_inheritance_for_report(
    v: Dict[str, Any],
    panel_mapper: Any,
    disease_name_to_inh: Dict[str, str],
) -> str:
    """Variant fields → gene-based panel → disease-name table from disease_gene_mapping.json."""
    inh = _coalesce_inheritance_for_report_variant(v)
    if inh:
        return inh
    inh = _inheritance_from_gene_panel_mapper(v.get("gene", ""), panel_mapper)
    if inh:
        return inh
    lbl = _disorder_label_for_template(v)
    inh = _inheritance_from_panel_disease_name_map(lbl, disease_name_to_inh)
    return inh


def _sanitize_diseases_for_report(
    v: Dict[str, Any], resolved_inheritance: str = ""
) -> List[Dict[str, Any]]:
    """Clean disease names in JSON; fills missing inheritance from resolved panel value."""
    out: List[Dict[str, Any]] = []
    for d in v.get("diseases") or []:
        if not isinstance(d, dict):
            continue
        raw = (d.get("name") or d.get("disease_name") or "").strip()
        if not raw:
            continue
        clean = _finalize_disorder_label(raw)
        if not clean:
            continue
        nd = dict(d)
        nd["name"] = clean
        nd["disease_name"] = clean
        if not (nd.get("inheritance") or "").strip() and resolved_inheritance:
            nd["inheritance"] = resolved_inheritance
        out.append(nd)
    return out


def _inheritance_display_for_report(v: Dict[str, Any]) -> str:
    """
    English phrase for PDF summary table (CN template maps these strings to Chinese).
    Codes from panel JSON: AR, AD, XL, etc.
    """
    code = _coalesce_inheritance_for_report_variant(v)
    if not code or code == "—":
        return "—"

    c0 = code.strip()
    cl = c0.lower()

    # Already a full phrase (e.g. localized text, or "Autosomal recessive")
    if len(c0) > 5 and (" " in c0 or not c0.isascii()):
        return c0

    mapping = {
        "ar": "Autosomal recessive",
        "ad": "Autosomal dominant",
        "xl": "X-linked",
        "x-linked": "X-linked",
        "xlinked": "X-linked",
        "xr": "X-linked recessive",
        "xd": "X-linked dominant",
        "mt": "Mitochondrial",
        "mitochondrial": "Mitochondrial",
        "pd": "Pseudoautosomal dominant",
        "pr": "Pseudoautosomal recessive",
        "ic": "Isolated cases",
        "mu": "Unknown",
        "somatic": "Somatic",
        "ol": "Oligogenic",
    }
    if cl in mapping:
        return mapping[cl]
    u = c0.upper()
    if u in ("AR", "AD", "XL", "XR", "XD", "MT"):
        return mapping.get(u.lower(), c0)
    return c0


def _disorder_label_for_template(v: Dict[str, Any]) -> str:
    diseases = v.get("diseases") or []
    if isinstance(diseases, list):
        for d in diseases:
            if isinstance(d, dict):
                name = (d.get("name") or d.get("disease_name") or "").strip()
                if name:
                    return _finalize_disorder_label(name)
            elif isinstance(d, str) and d.strip():
                return _finalize_disorder_label(d.strip())
    cd = (v.get("clinvar_disease") or v.get("clinvar_dn") or "").strip()
    if cd and cd not in (".",):
        one = _best_disorder_from_clinvar_field(cd)
        if one:
            return one
    hd = (v.get("hgmd_disease") or "").strip()
    if hd:
        one = _best_disorder_from_clinvar_field(hd) if ("|" in hd or ";" in hd) else _normalize_disorder_display(hd)
        if one:
            return one
    return "Unknown disorder"


def _variant_to_template_finding(v: Dict[str, Any]) -> Dict[str, Any]:
    """Shape expected by carrier_*.html / carrier_couples_*.html (orig/carrier_data_positive.json)."""
    disorder = _disorder_label_for_template(v)
    mutation = (
        (v.get("hgvsc") or "").strip()
        or (v.get("hgvsp") or "").strip()
        or (v.get("variant_id") or "")
    )
    rc = (v.get("reviewer_comment") or "").strip()
    parts = []
    if v.get("hgvsp"):
        parts.append(f"Protein: {v['hgvsp']}")
    if v.get("clinvar_sig"):
        parts.append(f"ClinVar: {v['clinvar_sig']}")
    auto_summary = rc if rc else (" ".join(parts) if parts else "")
    auto_gene_desc = ""
    if disorder and disorder != "Unknown disorder":
        auto_gene_desc = f"The {v.get('gene') or 'gene'} gene is associated with {disorder}."
    # Portal Generate Report tab: reviewer-edited text wins
    og = (v.get("report_gene_description") or "").strip()
    osum = (v.get("report_variant_summary") or "").strip()
    return {
        "gene": v.get("gene") or "",
        "mutation": mutation,
        "disorder": disorder,
        "inheritance": _inheritance_display_for_report(v),
        "classification": v.get("classification") or "VUS",
        "gene_description": og if og else auto_gene_desc,
        "variant_summary": osum if osum else auto_summary,
    }


def _finding_subject_raw(v: Dict[str, Any]) -> str:
    """Optional pipeline/portal keys to assign a variant to partner vs primary in couples mode."""
    s = (
        v.get("report_subject")
        or v.get("finding_subject")
        or v.get("couple_subject")
        or ""
    )
    s = str(s).strip().lower()
    if s in ("partner", "patient2", "p2", "spouse"):
        return "partner"
    return "primary"


def _split_template_findings(
    final_variants: List[Dict[str, Any]], is_couple: bool
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    if not is_couple:
        return ([_variant_to_template_finding(v) for v in final_variants], [])
    primary: List[Dict[str, Any]] = []
    partner: List[Dict[str, Any]] = []
    for v in final_variants:
        entry = _variant_to_template_finding(v)
        if _finding_subject_raw(v) == "partner":
            partner.append(entry)
        else:
            primary.append(entry)
    return primary, partner


# ══════════════════════════════════════════════════════════════
# Report JSON Generation
# ══════════════════════════════════════════════════════════════

def generate_report_json(
    order_id: str,
    sample_name: str,
    confirmed_variants: List[Dict[str, Any]],
    reviewer_info: Dict[str, Any],
    qc_summary: Dict[str, Any],
    output_dir: str,
    patient_info: Optional[Dict[str, Any]] = None,
    partner_info: Optional[Dict[str, Any]] = None,
    order_params: Optional[Dict[str, Any]] = None,
    report_language: str = "EN",
    disease_gene_json: Optional[str] = None,
    gene_knowledge_db: Optional[str] = None,
    gene_knowledge_enrich_on_report: bool = True,
    gene_knowledge_gemini_on_report: bool = True,
    gemini_api_key: Optional[str] = None,
    gene_knowledge_gemini_model: str = "gemini-2.5-flash",
    extra_result_json_paths: Optional[Sequence[str]] = None,
) -> str:
    """
    리뷰어 확정 결과를 바탕으로 report.json을 생성합니다.

    Args:
        order_id: 주문 ID
        sample_name: 샘플명
        confirmed_variants: 리뷰어가 확정한 변이 목록
            각 항목에는 reviewer_classification, reviewer_comment, reviewer_confirmed 포함
        reviewer_info: 리뷰어 정보 {name, id, institution, ...}
        qc_summary: QC 메트릭스 요약
        output_dir: 출력 디렉토리
        patient_info: 환자 정보 (선택)
        partner_info: 파트너 정보 (couple 검사 시, 선택)

    Returns:
        생성된 report.json 파일 경로

    disease_gene_json:
        Path to disease_gene_mapping.json — used to fill inheritance when missing on variants.

    gene_knowledge_db / Gemini:
        Optional SQLite + Gemini merge for **confirmed variants only** (review → Generate Report).
    """
    os.makedirs(output_dir, exist_ok=True)

    gk_db = (gene_knowledge_db or "").strip()
    if gk_db and gene_knowledge_enrich_on_report:
        from .gene_knowledge_db import enrich_confirmed_variants_for_report

        confirmed_variants = enrich_confirmed_variants_for_report(
            confirmed_variants,
            gene_knowledge_db=gk_db,
            gemini_api_key=(gemini_api_key or "").strip(),
            model=gene_knowledge_gemini_model,
            allow_gemini=bool(
                gene_knowledge_gemini_on_report and (gemini_api_key or "").strip()
            ),
        )

    panel_mapper: Any = None
    disease_name_to_inh: Dict[str, str] = {}
    if disease_gene_json and os.path.isfile(disease_gene_json):
        from .annotator import DiseaseGeneMapper

        panel_mapper = DiseaseGeneMapper(disease_gene_json)
        disease_name_to_inh = _load_panel_disease_inheritance_by_name(disease_gene_json)

    # 확정된 변이만 필터링
    final_variants = []
    for v in confirmed_variants:
        if not v.get("reviewer_confirmed", False):
            continue

        _clin = (v.get("clinvar_dn") or v.get("clinvar_disease") or "").strip()
        _inh = _resolve_inheritance_for_report(v, panel_mapper, disease_name_to_inh)
        _diseases = _sanitize_diseases_for_report(v, _inh)

        final_variant = {
            "variant_id": v.get("variant_id", ""),
            "chrom": v.get("chrom", ""),
            "pos": v.get("pos", ""),
            "ref": v.get("ref", ""),
            "alt": v.get("alt", ""),
            "gene": v.get("gene", ""),
            "transcript": v.get("transcript", ""),
            "clinical_nm": v.get("clinical_nm", ""),
            "hgvsc": v.get("hgvsc", ""),
            "hgvsp": v.get("hgvsp", ""),
            "effect": v.get("effect", ""),
            "zygosity": v.get("zygosity", ""),

            # 최종 분류 (리뷰어 확정)
            "classification": v.get("reviewer_classification") or v.get("acmg_classification", "VUS"),
            "acmg_criteria": v.get("acmg_criteria", []),

            # gnomAD
            "gnomad_af": v.get("gnomad_af"),

            # ClinVar (portal may send clinvar_dn or clinvar_disease)
            "clinvar_sig": v.get("clinvar_sig_primary", ""),
            "clinvar_dn": _clin,
            "clinvar_disease": _clin,

            # Disease-Gene Mapping (sanitized labels + coalesced inheritance for PDF / JSON)
            "diseases": _diseases,
            "inheritance": _inh,

            # HGMD
            "hgmd_class": v.get("hgmd_class", ""),
            "hgmd_disease": v.get("hgmd_disease", ""),

            # 리뷰어 코멘트
            "reviewer_comment": v.get("reviewer_comment", ""),

            # Couples: optional subject tag for partner vs primary (see _split_template_findings)
            "report_subject": v.get("report_subject")
            or v.get("finding_subject")
            or v.get("couple_subject"),

            # Portal Generate Report tab (editable before PDF)
            "report_gene_description": (v.get("report_gene_description") or "").strip(),
            "report_variant_summary": (v.get("report_variant_summary") or "").strip(),
        }
        final_variants.append(final_variant)

    # 질환별 그룹핑
    disease_groups = _group_by_disease(final_variants)

    # 캐리어 상태 판정
    carrier_status = _determine_carrier_status(final_variants)

    order_flat = _carrier_order_flat(order_params or {})

    # 환자 정보 기본값
    if patient_info is None:
        patient_info = {"name": sample_name}

    is_couple = bool(
        partner_info and (partner_info.get("name") or partner_info.get("dob"))
    )

    primary_findings, partner_findings = _split_template_findings(
        final_variants, is_couple
    )

    primary_patient = dict(patient_info)
    primary_patient.setdefault("name", sample_name)
    primary_patient["findings"] = primary_findings
    primary_patient.setdefault(
        "sample_id",
        order_flat.get("sample_id") or order_flat.get("medical_record_id") or "",
    )
    primary_patient.setdefault(
        "collection_date",
        order_flat.get("collection_date") or order_flat.get("sample_collection_date") or "",
    )

    partner_out: Optional[Dict[str, Any]] = None
    if partner_info:
        partner_out = dict(partner_info)
        partner_out["findings"] = partner_findings
        partner_out.setdefault(
            "sample_id",
            order_flat.get("partner_sample_id")
            or order_flat.get("patient2_sample_id")
            or "",
        )
        partner_out.setdefault(
            "collection_date",
            order_flat.get("partner_collection_date")
            or order_flat.get("patient2_collection_date")
            or primary_patient.get("collection_date")
            or "",
        )

    genes_n = order_flat.get("genes_evaluated_count")
    try:
        genes_evaluated_count = int(genes_n) if genes_n is not None and str(genes_n).strip() != "" else 302
    except (TypeError, ValueError):
        genes_evaluated_count = 302

    report_metadata = {
        "order_id": order_id,
        "report_date": _report_date_human_kst(),
        "pipeline_version": "gx-exome-v2.0",
        "report_version": "2.0",
        "is_couple": is_couple,
        "report_type": "Carrier Screening Report",
        "language": (report_language or "EN").strip().upper() or "EN",
        "hospital": order_flat.get("hospital_name")
        or order_flat.get("hospital")
        or "",
        "doctor": order_flat.get("doctor") or "",
    }

    report: Dict[str, Any] = {
        "version": "2.0",
        "type": "carrier_screening_report",
        "generated_at": now_kst_iso(),

        "report_metadata": report_metadata,

        # Jinja templates (carrier_*.html) read variants from primary_patient.findings
        "primary_patient": primary_patient,

        "partner": partner_out,

        "genes_evaluated_count": genes_evaluated_count,

        # 캐리어 상태 요약
        "carrier_status": carrier_status,

        # 확정된 변이 목록 (API / 프로그래밍용 — 템플릿은 findings 사용)
        "confirmed_variants": final_variants,
        "total_confirmed": len(final_variants),

        # 질환별 그룹
        "disease_groups": disease_groups,

        # QC 요약
        "qc_summary": {
            "mean_coverage": qc_summary.get("coverage", {}).get("mean_coverage"),
            "mapping_rate": qc_summary.get("alignment", {}).get("mapping_rate"),
            "total_reads": qc_summary.get("alignment", {}).get("total_reads"),
            "properly_paired_rate": qc_summary.get("alignment", {}).get("properly_paired_rate"),
        },

        # 리뷰어 정보
        "reviewer": reviewer_info,
    }

    # carrier_couples_*.html "Detailed Interpretations" uses data.findings (combined)
    if is_couple:
        report["findings"] = primary_findings + partner_findings

    # Dark genes supplement: read disk result.json → dark_genes_for_pdf → report.json only.
    # (Not copied from DB defaults; PDF uses detailed_sections + section_reviews when aligned.)
    # Same path fallbacks as generate_report_pdf when report output root ≠ pipeline output_dir.
    result_json_path = os.path.join(output_dir, "result.json")
    candidates = _dark_genes_result_json_paths_for_pdf(
        result_json_path, extra_result_json_paths
    )
    pdf_dg, src, any_file = _build_dark_genes_for_pdf_from_result_json_candidates(
        candidates, "generate_report_json"
    )
    if pdf_dg and src:
        report["dark_genes"] = pdf_dg
    elif any_file:
        logger.warning(
            "[generate_report_json] no usable dark_genes from candidates (newest-first) %s — "
            "supplement will be absent in PDF (check detailed_sections / approvals)",
            candidates,
        )

    output_path = os.path.join(output_dir, "report.json")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2, default=str)

    logger.info(
        f"Generated report.json: {output_path} "
        f"({len(final_variants)} confirmed variants, "
        f"carrier_status={carrier_status['status']})"
    )
    return output_path


def _group_by_disease(variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """변이를 질환별로 그룹핑합니다."""
    disease_map: Dict[str, List[Dict[str, Any]]] = {}

    for v in variants:
        # disease-gene mapping에서 질환 정보 추출
        diseases = v.get("diseases", [])
        clinvar_disease = v.get("clinvar_disease") or ""
        hgmd_disease = v.get("hgmd_disease") or ""

        disease_names = set()
        if diseases:
            for d in diseases:
                name = (d.get("name") or d.get("disease_name") or "") if isinstance(d, dict) else str(d)
                if name:
                    disease_names.add(name)
        if clinvar_disease and clinvar_disease not in (".", "not_provided", "not_specified"):
            disease_names.add(clinvar_disease)
        if hgmd_disease:
            disease_names.add(hgmd_disease)
        if not disease_names:
            disease_names.add("Unknown")

        for disease in disease_names:
            disease_map.setdefault(disease, []).append(v)

    groups = []
    for disease, vars_list in disease_map.items():
        pathogenic_count = sum(
            1 for v in vars_list
            if v.get("classification", "").lower() in ("pathogenic", "likely pathogenic")
        )
        groups.append({
            "disease": disease,
            "variants": vars_list,
            "variant_count": len(vars_list),
            "pathogenic_count": pathogenic_count,
        })

    # pathogenic 수가 많은 질환 순으로 정렬
    groups.sort(key=lambda g: g["pathogenic_count"], reverse=True)
    return groups


def _determine_carrier_status(variants: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    확정된 변이를 바탕으로 캐리어 상태를 판정합니다.

    판정 기준:
        - Pathogenic/Likely pathogenic 변이가 있으면 "carrier"
        - VUS만 있으면 "uncertain"
        - 없으면 "negative"
    """
    pathogenic_variants = [
        v for v in variants
        if v.get("classification", "").lower() in ("pathogenic", "likely pathogenic")
    ]
    vus_variants = [
        v for v in variants
        if v.get("classification", "").lower() == "vus"
    ]

    if pathogenic_variants:
        # 질환별 캐리어 상태
        carrier_diseases = set()
        for v in pathogenic_variants:
            diseases = v.get("diseases", [])
            clinvar_disease = v.get("clinvar_disease", "")
            if diseases:
                for d in diseases:
                    name = (d.get("name") or d.get("disease_name") or "") if isinstance(d, dict) else str(d)
                    if name:
                        carrier_diseases.add(name)
            elif clinvar_disease and clinvar_disease not in (".", "not_provided", "not_specified"):
                carrier_diseases.add(clinvar_disease)

        return {
            "status": "carrier",
            "pathogenic_count": len(pathogenic_variants),
            "carrier_diseases": list(carrier_diseases),
            "genes_affected": list(set(v.get("gene", "") for v in pathogenic_variants if v.get("gene"))),
            "summary": f"Carrier for {len(carrier_diseases)} disease(s) with {len(pathogenic_variants)} pathogenic variant(s)",
        }
    elif vus_variants:
        return {
            "status": "uncertain",
            "vus_count": len(vus_variants),
            "genes_affected": list(set(v.get("gene", "") for v in vus_variants if v.get("gene"))),
            "summary": f"{len(vus_variants)} variant(s) of uncertain significance found",
        }
    else:
        return {
            "status": "negative",
            "summary": "No pathogenic or likely pathogenic variants detected",
        }


# ══════════════════════════════════════════════════════════════
# PDF Report Generation (Jinja2 + WeasyPrint)
# Carrier_result/generate.py 패턴 통합
# ══════════════════════════════════════════════════════════════

def generate_report_pdf(
    report_json_path: str,
    output_dir: str,
    template_dir: Optional[str] = None,
    languages: Optional[List[str]] = None,
    extra_result_json_paths: Optional[Sequence[str]] = None,
) -> List[str]:
    """
    report.json을 바탕으로 다국어 PDF 리포트를 생성합니다.

    Carrier_result/generate.py의 Jinja2 + WeasyPrint 패턴을 따릅니다.
    template_dir에 Jinja2 HTML 템플릿이 있으면 사용하고,
    없으면 내장 HTML 템플릿을 사용합니다.

    Args:
        report_json_path: report.json 파일 경로
        output_dir: 출력 디렉토리
        template_dir: Jinja2 HTML 템플릿 디렉토리 (선택)
        languages: 생성할 언어 목록 (기본: ["EN"])
        extra_result_json_paths: 추가 ``result.json`` 후보 (파이프라인 출력 디렉토리 등)

    Returns:
        생성된 PDF 파일 경로 리스트
    """
    if languages is None:
        languages = ["EN"]

    try:
        with open(report_json_path, "r", encoding="utf-8") as f:
            report_data = json.load(f)
    except Exception as e:
        logger.error(f"Failed to read report JSON: {e}")
        return []

    _merge_dark_genes_from_result_json_for_pdf(
        report_data, report_json_path, extra_result_json_paths=extra_result_json_paths
    )

    try:
        from .dark_genes import sanitize_dark_genes_payload_for_pdf_render

        sanitize_dark_genes_payload_for_pdf_render(report_data)
    except Exception as e:
        logger.warning("[generate_report_pdf] dark_genes PDF sanitize skipped: %s", e)

    dg = report_data.get("dark_genes")
    if not isinstance(dg, dict) or not (
        (dg.get("report_detailed_html") or "").strip() or dg.get("status") == "error"
    ):
        logger.warning(
            "[generate_report_pdf] dark_genes missing or empty after merge/sanitize — "
            "PDF will have no supplementary hard-to-sequence block "
            "(check result.json dark_genes next to %s)",
            os.path.dirname(os.path.abspath(report_json_path)),
        )

    # Re-write report.json so on-disk JSON matches what WeasyPrint uses (merge is in-memory only above).
    try:
        with open(report_json_path, "w", encoding="utf-8") as f:
            json.dump(report_data, f, ensure_ascii=False, indent=2, default=str)
    except Exception as e:
        logger.warning("[generate_report_pdf] could not rewrite report.json after dark_genes merge: %s", e)

    os.makedirs(output_dir, exist_ok=True)

    if not template_dir:
        logger.error(
            "Carrier PDF: template_dir is None — output will use built-in HTML, not Jinja. "
            "Ensure data/carrier_report/carrier_EN.html exists or set template env vars."
        )

    is_couple = report_data.get("report_metadata", {}).get("is_couple", False)
    order_id = report_data.get("report_metadata", {}).get("order_id", "Unknown")
    raw_name = report_data.get("primary_patient", {}).get("name", "Patient")
    patient_name = "".join([c if c.isalnum() else "_" for c in raw_name]).replace("__", "_")

    generated_files = []

    for lang in languages:
        try:
            html_content = _render_html_for_language(
                report_data, lang, template_dir, is_couple
            )

            # HTML 저장
            html_filename = f"Report_{order_id}_{patient_name}_{lang}.html"
            html_path = os.path.join(output_dir, html_filename)
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_content)

            # PDF 생성
            pdf_filename = f"Report_{order_id}_{patient_name}_{lang}.pdf"
            pdf_path = os.path.join(output_dir, pdf_filename)

            try:
                from weasyprint import HTML as WeasyprintHTML
                base_url = template_dir if template_dir else output_dir
                WeasyprintHTML(string=html_content, base_url=base_url).write_pdf(pdf_path)
                generated_files.append(pdf_path)
                logger.info(f"Generated {lang} PDF: {pdf_filename}")
            except ImportError:
                logger.warning("weasyprint not installed, skipping PDF generation")
            except Exception as e:
                logger.error(f"Failed to generate PDF for {lang}: {e}")

        except Exception as e:
            logger.error(f"Error processing language {lang}: {e}")

    return generated_files


def _render_html_for_language(
    report_data: Dict[str, Any],
    lang: str,
    template_dir: Optional[str],
    is_couple: bool,
) -> str:
    """
    지정된 언어로 HTML 리포트를 렌더링합니다.

    template_dir에 Jinja2 템플릿이 있으면 사용하고,
    없으면 내장 HTML 생성 로직을 사용합니다.
    """
    lang_u = (lang or "EN").strip().upper() or "EN"

    if is_couple:
        template_name = f"carrier_couples_{lang_u}.html"
    else:
        template_name = f"carrier_{lang_u}.html"

    # Jinja2 템플릿 사용 시도
    if template_dir and os.path.isdir(template_dir):
        try:
            from jinja2 import Environment, FileSystemLoader

            env = Environment(loader=FileSystemLoader(template_dir))

            template = env.get_template(template_name)
            out = template.render(data=report_data)
            logger.info(
                "Carrier PDF: rendered Jinja template %s from %s",
                template_name,
                template_dir,
            )
            return out

        except Exception as e:
            logger.error(
                "Jinja2 carrier template failed (%s in %s): %s — using built-in HTML instead",
                template_name,
                template_dir,
                e,
                exc_info=True,
            )

    elif template_dir:
        logger.error(
            "Carrier PDF: template_dir is not a directory: %r — using built-in HTML",
            template_dir,
        )

    # 내장 HTML 생성
    return _render_builtin_html(report_data, lang_u)


def _render_builtin_html(report_data: Dict[str, Any], lang: str) -> str:
    """내장 HTML 템플릿으로 리포트를 렌더링합니다."""

    # 다국어 레이블
    labels = _get_labels(lang)

    sample_name = report_data.get("primary_patient", {}).get("name", "Unknown")
    order_id = report_data.get("report_metadata", {}).get("order_id", "Unknown")
    report_date = report_data.get("report_metadata", {}).get("report_date", "")
    carrier_status = report_data.get("carrier_status", {})
    confirmed_variants = report_data.get("confirmed_variants", [])
    qc_summary = report_data.get("qc_summary", {})
    disease_groups = report_data.get("disease_groups", [])
    partner = report_data.get("partner")

    # 캐리어 상태에 따른 색상
    status = carrier_status.get("status", "negative")
    status_color = {
        "carrier": "#dc3545",
        "uncertain": "#ffc107",
        "negative": "#28a745",
    }.get(status, "#6c757d")

    status_text = {
        "carrier": labels["status_carrier"],
        "uncertain": labels["status_uncertain"],
        "negative": labels["status_negative"],
    }.get(status, status)

    # 변이 테이블 행 생성
    variant_rows = ""
    for v in confirmed_variants:
        classification = v.get("classification", "VUS")
        cls_color = _classification_color(classification)

        gnomad_af = v.get("gnomad_af")
        af_str = f"{gnomad_af:.6f}" if gnomad_af is not None else "N/A"

        # 질환 정보
        disease_str = ""
        diseases = v.get("diseases", [])
        if diseases:
            disease_str = "; ".join(
                (d.get("name") or d.get("disease_name") or "")
                for d in diseases
                if isinstance(d, dict) and (d.get("name") or d.get("disease_name"))
            )
        if not disease_str:
            disease_str = v.get("clinvar_disease", "")

        variant_rows += f"""
        <tr>
            <td>{v.get('gene', '')}</td>
            <td>{v.get('chrom', '')}:{v.get('pos', '')}</td>
            <td>{v.get('hgvsc', '')}</td>
            <td>{v.get('hgvsp', '')}</td>
            <td>{v.get('zygosity', '')}</td>
            <td style="color: {cls_color}; font-weight: bold;">{classification}</td>
            <td>{af_str}</td>
            <td>{v.get('inheritance', '')}</td>
            <td>{disease_str}</td>
            <td>{v.get('reviewer_comment', '')}</td>
        </tr>
        """

    # 질환 요약 섹션
    disease_summary_html = ""
    if disease_groups:
        disease_summary_html = f"<h2>{labels['disease_summary']}</h2><ul>"
        for dg in disease_groups:
            disease_summary_html += (
                f"<li><strong>{dg['disease']}</strong>: "
                f"{dg['variant_count']} variant(s), "
                f"{dg['pathogenic_count']} pathogenic</li>"
            )
        disease_summary_html += "</ul>"

    # 파트너 정보 섹션
    partner_section = ""
    if partner:
        partner_section = f"""
        <div class="partner-info">
            <h3>{labels['partner_info']}</h3>
            <p><strong>{labels['name']}:</strong> {partner.get('name', 'N/A')}</p>
        </div>
        """

    html = f"""<!DOCTYPE html>
<html lang="{lang.lower()}">
<head>
    <meta charset="UTF-8">
    <title>{labels['title']} - {sample_name}</title>
    <style>
        @page {{
            size: A4;
            margin: 2cm;
        }}
        body {{
            font-family: 'Noto Sans KR', 'Noto Sans SC', Arial, sans-serif;
            margin: 0;
            padding: 40px;
            color: #333;
            font-size: 10pt;
            line-height: 1.5;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
            font-size: 18pt;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            font-size: 14pt;
            border-bottom: 1px solid #bdc3c7;
            padding-bottom: 5px;
        }}
        .header-info {{
            display: flex;
            justify-content: space-between;
            margin-bottom: 20px;
        }}
        .header-info div {{
            flex: 1;
        }}
        .status-badge {{
            display: inline-block;
            padding: 8px 20px;
            border-radius: 20px;
            color: white;
            font-weight: bold;
            font-size: 1.1em;
            background-color: {status_color};
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            font-size: 0.85em;
        }}
        th {{
            background-color: #3498db;
            color: white;
            padding: 8px 6px;
            text-align: left;
            font-size: 0.85em;
        }}
        td {{
            padding: 6px;
            border-bottom: 1px solid #ecf0f1;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        .qc-grid {{
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 15px;
            margin: 15px 0;
        }}
        .qc-item {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            border: 1px solid #e9ecef;
        }}
        .qc-item .value {{
            font-size: 1.3em;
            font-weight: bold;
            color: #2c3e50;
        }}
        .qc-item .label {{
            font-size: 0.8em;
            color: #7f8c8d;
            margin-top: 5px;
        }}
        .partner-info {{
            background: #f0f8ff;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            font-size: 0.8em;
            color: #7f8c8d;
        }}
        .disclaimer {{
            background: #fff3cd;
            padding: 12px;
            border-radius: 6px;
            margin: 20px 0;
            font-size: 0.85em;
        }}
        @media print {{
            body {{ margin: 0; padding: 20px; }}
            .no-print {{ display: none; }}
        }}
    </style>
</head>
<body>
    <h1>{labels['title']}</h1>

    <div class="header-info">
        <div>
            <p><strong>{labels['sample']}:</strong> {sample_name}</p>
            <p><strong>{labels['order_id']}:</strong> {order_id}</p>
            <p><strong>{labels['report_date']}:</strong> {report_date}</p>
        </div>
        <div style="text-align: right;">
            <p><strong>{labels['result']}:</strong></p>
            <span class="status-badge">{status_text}</span>
        </div>
    </div>

    {partner_section}

    <h2>{labels['qc_summary']}</h2>
    <div class="qc-grid">
        <div class="qc-item">
            <div class="value">{_fmt_val(qc_summary.get('mean_coverage'), '.1f')}x</div>
            <div class="label">{labels['mean_coverage']}</div>
        </div>
        <div class="qc-item">
            <div class="value">{_fmt_val(qc_summary.get('mapping_rate'), '.1f')}%</div>
            <div class="label">{labels['mapping_rate']}</div>
        </div>
        <div class="qc-item">
            <div class="value">{_fmt_reads(qc_summary.get('total_reads'))}</div>
            <div class="label">{labels['total_reads']}</div>
        </div>
        <div class="qc-item">
            <div class="value">{_fmt_val(qc_summary.get('properly_paired_rate'), '.1f')}%</div>
            <div class="label">{labels['properly_paired']}</div>
        </div>
    </div>

    {disease_summary_html}

    <h2>{labels['confirmed_variants']} ({len(confirmed_variants)})</h2>
    <table>
        <thead>
            <tr>
                <th>{labels['gene']}</th>
                <th>{labels['position']}</th>
                <th>HGVS.c</th>
                <th>HGVS.p</th>
                <th>{labels['zygosity']}</th>
                <th>{labels['classification']}</th>
                <th>gnomAD AF</th>
                <th>{labels['inheritance']}</th>
                <th>{labels['disease']}</th>
                <th>{labels['comment']}</th>
            </tr>
        </thead>
        <tbody>
            {variant_rows if variant_rows else f'<tr><td colspan="10" style="text-align:center;">{labels["no_variants"]}</td></tr>'}
        </tbody>
    </table>

    {_builtin_dark_genes_supplement_html(report_data)}

    <div class="disclaimer">
        {labels['disclaimer']}
    </div>

    <div class="footer">
        <p>{labels['footer_pipeline']}</p>
        <p>{labels['footer_acmg']}</p>
        <p>{labels['reviewer']}: {(report_data.get('reviewer') or {}).get('name', 'N/A')}</p>
    </div>
</body>
</html>"""

    return html


def _builtin_dark_genes_supplement_html(report_data: Dict[str, Any]) -> str:
    """Built-in carrier HTML (no Jinja): same supplement as carrier_EN.html when HTML is present."""
    dg = report_data.get("dark_genes")
    if not isinstance(dg, dict):
        return ""
    if dg.get("status") == "error" and (dg.get("report_summary") or "").strip():
        return (
            '<h2 style="page-break-before: always;">Supplementary analysis (hard-to-sequence regions)</h2>'
            f'<p style="white-space:pre-wrap;font-size:9pt;">'
            f"{html.escape(str(dg.get('report_summary') or ''), quote=True)}</p>"
        )
    html_snip = (dg.get("report_detailed_html") or "").strip()
    if not html_snip:
        return ""
    return (
        '<h2 style="page-break-before: always;">Supplementary analysis (hard-to-sequence regions)</h2>'
        '<div style="font-size:9pt;color:#334155;">' + html_snip + "</div>"
    )


def _get_labels(lang: str) -> Dict[str, str]:
    """다국어 레이블을 반환합니다."""
    labels = {
        "EN": {
            "title": "Carrier Screening Report",
            "sample": "Sample",
            "order_id": "Order ID",
            "report_date": "Report Date",
            "result": "Result",
            "status_carrier": "Carrier Detected",
            "status_uncertain": "Uncertain Significance",
            "status_negative": "Negative",
            "qc_summary": "QC Summary",
            "mean_coverage": "Mean Coverage",
            "mapping_rate": "Mapping Rate",
            "total_reads": "Total Reads",
            "properly_paired": "Properly Paired",
            "disease_summary": "Disease Summary",
            "confirmed_variants": "Confirmed Variants",
            "gene": "Gene",
            "position": "Position",
            "zygosity": "Zygosity",
            "classification": "Classification",
            "inheritance": "Inheritance",
            "disease": "Disease",
            "comment": "Comment",
            "no_variants": "No confirmed variants",
            "partner_info": "Partner Information",
            "name": "Name",
            "disclaimer": (
                "This report is intended for use by qualified healthcare professionals. "
                "Variant classification follows ACMG/AMP guidelines. All results should be "
                "interpreted in the context of clinical findings and family history."
            ),
            "footer_pipeline": "Generated by Carrier Screening Pipeline v2.0",
            "footer_acmg": "Variant classification follows ACMG/AMP 2015 guidelines.",
            "reviewer": "Reviewed by",
        },
        "CN": {
            "title": "携带者筛查报告",
            "sample": "样本",
            "order_id": "订单号",
            "report_date": "报告日期",
            "result": "结果",
            "status_carrier": "检测到携带者",
            "status_uncertain": "意义不明",
            "status_negative": "阴性",
            "qc_summary": "质控摘要",
            "mean_coverage": "平均覆盖度",
            "mapping_rate": "比对率",
            "total_reads": "总读数",
            "properly_paired": "正确配对",
            "disease_summary": "疾病摘要",
            "confirmed_variants": "确认变异",
            "gene": "基因",
            "position": "位置",
            "zygosity": "合子性",
            "classification": "分类",
            "inheritance": "遗传方式",
            "disease": "疾病",
            "comment": "备注",
            "no_variants": "无确认变异",
            "partner_info": "伴侣信息",
            "name": "姓名",
            "disclaimer": (
                "本报告仅供合格医疗专业人员使用。"
                "变异分类遵循ACMG/AMP指南。所有结果应结合临床发现和家族史进行解读。"
            ),
            "footer_pipeline": "由携带者筛查管线 v2.0 生成",
            "footer_acmg": "变异分类遵循ACMG/AMP 2015指南。",
            "reviewer": "审核人",
        },
        "KO": {
            "title": "보인자 선별검사 보고서",
            "sample": "검체",
            "order_id": "주문번호",
            "report_date": "보고일",
            "result": "결과",
            "status_carrier": "보인자 검출",
            "status_uncertain": "임상적 의의 불확실",
            "status_negative": "음성",
            "qc_summary": "QC 요약",
            "mean_coverage": "평균 커버리지",
            "mapping_rate": "매핑률",
            "total_reads": "총 리드 수",
            "properly_paired": "정상 페어링",
            "disease_summary": "질환 요약",
            "confirmed_variants": "확정 변이",
            "gene": "유전자",
            "position": "위치",
            "zygosity": "접합성",
            "classification": "분류",
            "inheritance": "유전양식",
            "disease": "질환",
            "comment": "코멘트",
            "no_variants": "확정된 변이 없음",
            "partner_info": "파트너 정보",
            "name": "이름",
            "disclaimer": (
                "본 보고서는 자격을 갖춘 의료 전문가를 위한 것입니다. "
                "변이 분류는 ACMG/AMP 가이드라인을 따릅니다. 모든 결과는 "
                "임상 소견 및 가족력을 고려하여 해석되어야 합니다."
            ),
            "footer_pipeline": "Carrier Screening Pipeline v2.0으로 생성됨",
            "footer_acmg": "변이 분류는 ACMG/AMP 2015 가이드라인을 따릅니다.",
            "reviewer": "검토자",
        },
    }

    return labels.get(lang.upper(), labels["EN"])


def _classification_color(classification: str) -> str:
    """ACMG 분류에 따른 색상을 반환합니다."""
    return {
        "Pathogenic": "#dc3545",
        "Likely pathogenic": "#fd7e14",
        "VUS": "#ffc107",
        "Likely benign": "#20c997",
        "Benign": "#28a745",
    }.get(classification, "#6c757d")


def _fmt_val(val, fmt: str = "") -> str:
    """값을 포맷팅합니다."""
    if val is None:
        return "N/A"
    try:
        return f"{val:{fmt}}" if fmt else str(val)
    except (ValueError, TypeError):
        return str(val)


def _fmt_reads(val) -> str:
    """리드 수를 포맷팅합니다 (예: 12.3M)."""
    if val is None:
        return "N/A"
    try:
        n = int(val)
        if n >= 1_000_000:
            return f"{n / 1_000_000:.1f}M"
        elif n >= 1_000:
            return f"{n / 1_000:.1f}K"
        return str(n)
    except (ValueError, TypeError):
        return str(val)
