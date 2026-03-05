"""
Report Generator Module

리뷰어가 확정한 결과를 바탕으로 최종 리포트를 생성합니다.

생성 파일:
    - report.json       : 최종 확정 결과 (Portal → 고객 전달용)
    - report.pdf        : PDF 리포트 (선택적)

report.json 구조:
    - 리뷰어가 확정한 변이 목록
    - 각 변이의 ACMG 분류 (리뷰어 확정 버전)
    - 리뷰어 코멘트
    - 질환 관련 정보
    - QC 요약
"""

import os
import json
import logging
from datetime import datetime
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)


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

            # 리뷰어 코멘트
            "reviewer_comment": v.get("reviewer_comment", ""),
        }
        final_variants.append(final_variant)

    # 질환별 그룹핑
    disease_groups = _group_by_disease(final_variants)

    # 캐리어 상태 판정
    carrier_status = _determine_carrier_status(final_variants)

    report = {
        "version": "1.0",
        "type": "carrier_screening_report",
        "generated_at": datetime.now().isoformat(),
        "order_id": order_id,
        "sample_name": sample_name,

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
        },

        # 리뷰어 정보
        "reviewer": reviewer_info,

        # 메타데이터
        "metadata": {
            "pipeline_version": "carrier-screening-v1.0",
            "report_version": "1.0",
        },
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
        disease = v.get("clinvar_disease") or "Unknown"
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
            disease = v.get("clinvar_disease") or "Unknown"
            if disease != "Unknown":
                carrier_diseases.add(disease)

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
# PDF Report Generation (Optional)
# ══════════════════════════════════════════════════════════════

def generate_report_pdf(
    report_json_path: str,
    output_dir: str,
    template_dir: Optional[str] = None,
) -> Optional[str]:
    """
    report.json을 바탕으로 PDF 리포트를 생성합니다.

    Args:
        report_json_path: report.json 파일 경로
        output_dir: 출력 디렉토리
        template_dir: HTML 템플릿 디렉토리 (선택)

    Returns:
        생성된 PDF 파일 경로 또는 None
    """
    try:
        with open(report_json_path, "r", encoding="utf-8") as f:
            report_data = json.load(f)
    except Exception as e:
        logger.error(f"Failed to read report JSON: {e}")
        return None

    os.makedirs(output_dir, exist_ok=True)

    try:
        # HTML 리포트 생성
        html_content = _render_report_html(report_data, template_dir)
        html_path = os.path.join(output_dir, "report.html")
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        # HTML → PDF 변환
        pdf_path = os.path.join(output_dir, "report.pdf")
        try:
            from weasyprint import HTML
            HTML(string=html_content).write_pdf(pdf_path)
            logger.info(f"Generated report PDF: {pdf_path}")
            return pdf_path
        except ImportError:
            logger.warning("weasyprint not installed, skipping PDF generation")
            return None

    except Exception as e:
        logger.error(f"Failed to generate PDF report: {e}")
        return None


def _render_report_html(report_data: Dict[str, Any], template_dir: Optional[str] = None) -> str:
    """report.json 데이터를 HTML로 렌더링합니다."""
    sample_name = report_data.get("sample_name", "Unknown")
    order_id = report_data.get("order_id", "Unknown")
    generated_at = report_data.get("generated_at", "")
    carrier_status = report_data.get("carrier_status", {})
    confirmed_variants = report_data.get("confirmed_variants", [])
    qc_summary = report_data.get("qc_summary", {})

    # 캐리어 상태에 따른 색상
    status = carrier_status.get("status", "negative")
    status_color = {
        "carrier": "#dc3545",
        "uncertain": "#ffc107",
        "negative": "#28a745",
    }.get(status, "#6c757d")

    # 변이 테이블 행 생성
    variant_rows = ""
    for v in confirmed_variants:
        classification = v.get("classification", "VUS")
        cls_color = {
            "Pathogenic": "#dc3545",
            "Likely pathogenic": "#fd7e14",
            "VUS": "#ffc107",
            "Likely benign": "#20c997",
            "Benign": "#28a745",
        }.get(classification, "#6c757d")

        gnomad_af = v.get("gnomad_af")
        af_str = f"{gnomad_af:.6f}" if gnomad_af is not None else "N/A"

        variant_rows += f"""
        <tr>
            <td>{v.get('gene', '')}</td>
            <td>{v.get('chrom', '')}:{v.get('pos', '')}</td>
            <td>{v.get('hgvsc', '')}</td>
            <td>{v.get('hgvsp', '')}</td>
            <td>{v.get('zygosity', '')}</td>
            <td style="color: {cls_color}; font-weight: bold;">{classification}</td>
            <td>{af_str}</td>
            <td>{v.get('clinvar_disease', '')}</td>
            <td>{v.get('reviewer_comment', '')}</td>
        </tr>
        """

    html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
    <meta charset="UTF-8">
    <title>Carrier Screening Report - {sample_name}</title>
    <style>
        body {{ font-family: 'Noto Sans KR', Arial, sans-serif; margin: 40px; color: #333; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .header-info {{ display: flex; justify-content: space-between; margin-bottom: 20px; }}
        .header-info div {{ flex: 1; }}
        .status-badge {{
            display: inline-block; padding: 8px 20px; border-radius: 20px;
            color: white; font-weight: bold; font-size: 1.1em;
            background-color: {status_color};
        }}
        table {{ width: 100%; border-collapse: collapse; margin: 15px 0; font-size: 0.9em; }}
        th {{ background-color: #3498db; color: white; padding: 10px 8px; text-align: left; }}
        td {{ padding: 8px; border-bottom: 1px solid #ddd; }}
        tr:nth-child(even) {{ background-color: #f8f9fa; }}
        .qc-grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; margin: 15px 0; }}
        .qc-item {{ background: #f8f9fa; padding: 15px; border-radius: 8px; text-align: center; }}
        .qc-item .value {{ font-size: 1.5em; font-weight: bold; color: #2c3e50; }}
        .qc-item .label {{ font-size: 0.85em; color: #7f8c8d; }}
        .footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd;
                   font-size: 0.85em; color: #7f8c8d; }}
        @media print {{ body {{ margin: 20px; }} }}
    </style>
</head>
<body>
    <h1>Carrier Screening Report</h1>

    <div class="header-info">
        <div>
            <p><strong>Sample:</strong> {sample_name}</p>
            <p><strong>Order ID:</strong> {order_id}</p>
            <p><strong>Report Date:</strong> {generated_at[:10] if generated_at else ''}</p>
        </div>
        <div style="text-align: right;">
            <p><strong>Result:</strong></p>
            <span class="status-badge">{carrier_status.get('summary', '')}</span>
        </div>
    </div>

    <h2>QC Summary</h2>
    <div class="qc-grid">
        <div class="qc-item">
            <div class="value">{qc_summary.get('mean_coverage', 'N/A')}</div>
            <div class="label">Mean Coverage</div>
        </div>
        <div class="qc-item">
            <div class="value">{qc_summary.get('mapping_rate', 'N/A')}%</div>
            <div class="label">Mapping Rate</div>
        </div>
        <div class="qc-item">
            <div class="value">{qc_summary.get('total_reads', 'N/A')}</div>
            <div class="label">Total Reads</div>
        </div>
    </div>

    <h2>Confirmed Variants ({len(confirmed_variants)})</h2>
    <table>
        <thead>
            <tr>
                <th>Gene</th>
                <th>Position</th>
                <th>HGVS.c</th>
                <th>HGVS.p</th>
                <th>Zygosity</th>
                <th>Classification</th>
                <th>gnomAD AF</th>
                <th>Disease</th>
                <th>Comment</th>
            </tr>
        </thead>
        <tbody>
            {variant_rows if variant_rows else '<tr><td colspan="9" style="text-align:center;">No confirmed variants</td></tr>'}
        </tbody>
    </table>

    <div class="footer">
        <p>This report was generated by the Carrier Screening Pipeline v1.0.</p>
        <p>Variant classification follows ACMG/AMP guidelines. All results should be interpreted
           by a qualified clinical geneticist in the context of clinical findings.</p>
        <p>Reviewer: {report_data.get('reviewer', {}).get('name', 'N/A')}</p>
    </div>
</body>
</html>"""

    return html
