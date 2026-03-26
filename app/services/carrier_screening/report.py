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
    - version, type, generated_at
    - order_id, sample_name
    - primary_patient: 환자 정보
    - partner: 파트너 정보 (couple 검사 시)
    - carrier_status: 캐리어 상태 판정
    - confirmed_variants: 리뷰어 확정 변이 목록
    - disease_groups: 질환별 변이 그룹
    - qc_summary: QC 메트릭스 요약
    - reviewer: 리뷰어 정보
    - report_metadata: 리포트 메타데이터
"""

import os
import json
import logging
from typing import Dict, Any, List, Optional

from ...datetime_kst import now_kst_date_iso, now_kst_iso

logger = logging.getLogger(__name__)


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
    """
    os.makedirs(output_dir, exist_ok=True)

    # 확정된 변이만 필터링
    final_variants = []
    for v in confirmed_variants:
        if not v.get("reviewer_confirmed", False):
            continue

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

            # ClinVar
            "clinvar_sig": v.get("clinvar_sig_primary", ""),
            "clinvar_disease": v.get("clinvar_dn", ""),

            # Disease-Gene Mapping
            "diseases": v.get("diseases", []),
            "inheritance": v.get("inheritance", ""),

            # HGMD
            "hgmd_class": v.get("hgmd_class", ""),
            "hgmd_disease": v.get("hgmd_disease", ""),

            # 리뷰어 코멘트
            "reviewer_comment": v.get("reviewer_comment", ""),
        }
        final_variants.append(final_variant)

    # 질환별 그룹핑
    disease_groups = _group_by_disease(final_variants)

    # 캐리어 상태 판정
    carrier_status = _determine_carrier_status(final_variants)

    # 환자 정보 기본값
    if patient_info is None:
        patient_info = {"name": sample_name}

    is_couple = bool(
        partner_info and (partner_info.get("name") or partner_info.get("dob"))
    )

    report = {
        "version": "2.0",
        "type": "carrier_screening_report",
        "generated_at": now_kst_iso(),

        # 리포트 메타데이터 (Carrier_result/generate.py 호환)
        "report_metadata": {
            "order_id": order_id,
            "report_date": now_kst_date_iso(),
            "pipeline_version": "carrier-screening-v2.0",
            "report_version": "2.0",
            "is_couple": is_couple,
        },

        # 환자 정보
        "primary_patient": patient_info,

        # 파트너 정보 (couple 검사 시)
        "partner": partner_info,

        # 캐리어 상태 요약
        "carrier_status": carrier_status,

        # 확정된 변이 목록
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
                name = d.get("name", "")
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
                    name = d.get("name", "")
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

    os.makedirs(output_dir, exist_ok=True)

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
    # Jinja2 템플릿 사용 시도
    if template_dir and os.path.isdir(template_dir):
        try:
            from jinja2 import Environment, FileSystemLoader

            env = Environment(loader=FileSystemLoader(template_dir))

            # Carrier_result/generate.py 패턴: carrier_{lang}.html / carrier_couples_{lang}.html
            if is_couple:
                template_name = f"carrier_couples_{lang}.html"
            else:
                template_name = f"carrier_{lang}.html"

            template = env.get_template(template_name)
            return template.render(data=report_data)

        except Exception as e:
            logger.warning(f"Jinja2 template rendering failed ({e}), falling back to built-in HTML")

    # 내장 HTML 생성
    return _render_builtin_html(report_data, lang)


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
            disease_str = "; ".join(d.get("name", "") for d in diseases if d.get("name"))
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
