"""
Review Data Generator Module

파이프라인 결과를 기반으로 Portal 리뷰 페이지에 필요한 데이터를 생성합니다.

생성 파일:
    - result.json       : 변이 목록 + annotation + ACMG 분류 (Portal 리뷰 페이지용)
    - qc_summary.json   : QC 메트릭스 요약 (coverage, mapping rate 등)
    - variants.tsv      : 변이 목록 TSV (백업/다운로드용)
    - IGV snapshots     : 주요 변이 위치의 IGV 스냅샷 이미지 (파이프라인에서 생성)
"""

import os
import re
import csv
import json
import glob
import logging
from datetime import datetime
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# QC Summary Extraction
# ══════════════════════════════════════════════════════════════

def extract_qc_summary(analysis_dir: str, sample_name: str) -> Dict[str, Any]:
    """
    파이프라인 분석 디렉토리에서 QC 메트릭스를 추출합니다.

    Nextflow 파이프라인이 생성하는 주요 QC 파일:
        - {sample}.mosdepth.summary.txt   (coverage)
        - {sample}.flagstat               (alignment stats)
        - {sample}.stats                  (samtools stats)
        - {sample}.coverage_report.txt    (target coverage)
    """
    qc = {
        "sample_name": sample_name,
        "generated_at": datetime.now().isoformat(),
        "coverage": {},
        "alignment": {},
        "variant_stats": {},
    }

    # mosdepth coverage
    mosdepth_files = glob.glob(os.path.join(analysis_dir, "**", f"*mosdepth*summary*"), recursive=True)
    for mf in mosdepth_files:
        try:
            qc["coverage"].update(_parse_mosdepth_summary(mf))
        except Exception as e:
            logger.warning(f"Failed to parse mosdepth file {mf}: {e}")

    # flagstat alignment
    flagstat_files = glob.glob(os.path.join(analysis_dir, "**", f"*flagstat*"), recursive=True)
    for ff in flagstat_files:
        try:
            qc["alignment"].update(_parse_flagstat(ff))
        except Exception as e:
            logger.warning(f"Failed to parse flagstat file {ff}: {e}")

    # samtools stats
    stats_files = glob.glob(os.path.join(analysis_dir, "**", f"*{sample_name}*.stats"), recursive=True)
    for sf in stats_files:
        try:
            parsed = _parse_samtools_stats(sf)
            if parsed:
                qc["alignment"].update(parsed)
        except Exception as e:
            logger.warning(f"Failed to parse stats file {sf}: {e}")

    # coverage report (target-specific)
    cov_report_files = glob.glob(os.path.join(analysis_dir, "**", f"*coverage_report*"), recursive=True)
    for cf in cov_report_files:
        try:
            qc["coverage"].update(_parse_coverage_report(cf))
        except Exception as e:
            logger.warning(f"Failed to parse coverage report {cf}: {e}")

    return qc


def _parse_mosdepth_summary(path: str) -> Dict[str, Any]:
    """mosdepth summary.txt 파싱"""
    result = {}
    with open(path, "r") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < len(header):
                continue
            row = dict(zip(header, parts))
            region = row.get("chrom", "")
            if region == "total" or region == "total_region":
                try:
                    result["mean_coverage"] = float(row.get("mean", 0))
                    result["min_coverage"] = float(row.get("min", 0))
                    result["max_coverage"] = float(row.get("max", 0))
                except (ValueError, TypeError):
                    pass
    return result


def _parse_flagstat(path: str) -> Dict[str, Any]:
    """samtools flagstat 파싱"""
    result = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if "mapped" in line.lower() and "%" in line:
                m = re.search(r"(\d+)\s+\+\s+\d+\s+mapped\s+\((\d+\.?\d*)%", line)
                if m:
                    result["mapped_reads"] = int(m.group(1))
                    result["mapping_rate"] = float(m.group(2))
            elif "total" in line.lower() and "QC-passed" in line:
                m = re.search(r"(\d+)\s+\+\s+\d+\s+in total", line)
                if m:
                    result["total_reads"] = int(m.group(1))
            elif "properly paired" in line.lower():
                m = re.search(r"(\d+)\s+\+\s+\d+\s+properly paired\s+\((\d+\.?\d*)%", line)
                if m:
                    result["properly_paired"] = int(m.group(1))
                    result["properly_paired_rate"] = float(m.group(2))
            elif "duplicates" in line.lower():
                m = re.search(r"(\d+)\s+\+\s+\d+\s+duplicates", line)
                if m:
                    result["duplicates"] = int(m.group(1))
    return result


def _parse_samtools_stats(path: str) -> Dict[str, Any]:
    """samtools stats 파싱"""
    result = {}
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("SN\t"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            key = parts[1].rstrip(":")
            val = parts[2]
            if key == "raw total sequences":
                result["raw_total_sequences"] = int(val)
            elif key == "reads mapped":
                result.setdefault("mapped_reads", int(val))
            elif key == "insert size average":
                result["insert_size_avg"] = float(val)
            elif key == "insert size standard deviation":
                result["insert_size_std"] = float(val)
            elif key == "average quality":
                result["average_quality"] = float(val)
    return result


def _parse_coverage_report(path: str) -> Dict[str, Any]:
    """target coverage report 파싱"""
    result = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if ":" in line:
                key, _, val = line.partition(":")
                key = key.strip().lower().replace(" ", "_")
                val = val.strip()
                try:
                    if "%" in val:
                        result[key] = float(val.replace("%", ""))
                    elif "." in val:
                        result[key] = float(val)
                    else:
                        result[key] = int(val)
                except (ValueError, TypeError):
                    result[key] = val
    return result


# ══════════════════════════════════════════════════════════════
# Result JSON Generation
# ══════════════════════════════════════════════════════════════

def generate_result_json(
    annotated_variants: List[Dict[str, Any]],
    acmg_results: List[Dict[str, Any]],
    qc_summary: Dict[str, Any],
    sample_name: str,
    order_id: str,
    disease_bed_info: Optional[Dict[str, Any]] = None,
    output_dir: str = "",
) -> str:
    """
    Portal 리뷰 페이지용 result.json을 생성합니다.

    Args:
        annotated_variants: annotator.annotate() 결과 리스트
        acmg_results: acmg.classify_variant() 결과 리스트
        qc_summary: extract_qc_summary() 결과
        sample_name: 샘플명
        order_id: 주문 ID
        disease_bed_info: 질환 BED 정보 (선택)
        output_dir: 출력 디렉토리

    Returns:
        생성된 result.json 파일 경로
    """
    os.makedirs(output_dir, exist_ok=True)

    # 변이 + ACMG 결합
    variants_for_review = []
    for i, (var, acmg) in enumerate(zip(annotated_variants, acmg_results)):
        entry = {
            "variant_id": f"VAR_{i+1:04d}",
            **var,
            "acmg_classification": acmg.get("final_classification", "VUS"),
            "acmg_criteria": acmg.get("final_criteria", []),
            "acmg_reasoning": acmg.get("final_reasoning", ""),
            "acmg_rule_based": acmg.get("rule_based", {}),
            "acmg_ai": acmg.get("ai"),
            # 리뷰어 입력 필드 (Portal에서 채움)
            "reviewer_classification": None,
            "reviewer_comment": None,
            "reviewer_confirmed": False,
        }
        variants_for_review.append(entry)

    # 변이 통계
    variant_stats = _compute_variant_stats(variants_for_review)

    # IGV 스냅샷 경로 매핑
    igv_snapshots = _find_igv_snapshots(output_dir, variants_for_review)

    result = {
        "version": "1.0",
        "type": "carrier_screening_result",
        "generated_at": datetime.now().isoformat(),
        "order_id": order_id,
        "sample_name": sample_name,
        "status": "pending_review",

        # QC 요약
        "qc_summary": qc_summary,

        # 변이 통계
        "variant_stats": variant_stats,

        # 변이 목록 (리뷰 대상)
        "variants": variants_for_review,

        # 질환 BED 정보
        "disease_panel": disease_bed_info or {},

        # IGV 스냅샷 매핑
        "igv_snapshots": igv_snapshots,

        # 메타데이터
        "metadata": {
            "pipeline_version": "carrier-screening-v1.0",
            "annotation_databases": {
                "clinvar": True,
                "gnomad": True,
                "snpeff": True,
                "clingen": True,
                "mane": True,
            },
            "acmg_mode": "rule_based",
        },
    }

    output_path = os.path.join(output_dir, "result.json")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2, default=str)

    logger.info(f"Generated result.json: {output_path} ({len(variants_for_review)} variants)")
    return output_path


def _compute_variant_stats(variants: List[Dict[str, Any]]) -> Dict[str, Any]:
    """변이 통계를 계산합니다."""
    total = len(variants)
    if total == 0:
        return {"total": 0}

    classifications = {}
    zygosities = {}
    effects = {}

    for v in variants:
        cls = v.get("acmg_classification", "VUS")
        classifications[cls] = classifications.get(cls, 0) + 1

        zyg = v.get("zygosity", "Unknown")
        zygosities[zyg] = zygosities.get(zyg, 0) + 1

        eff = v.get("effect", "Unknown")
        effects[eff] = effects.get(eff, 0) + 1

    pathogenic_count = sum(
        classifications.get(c, 0)
        for c in ("Pathogenic", "Likely pathogenic")
    )

    return {
        "total": total,
        "pathogenic_or_likely": pathogenic_count,
        "vus": classifications.get("VUS", 0),
        "benign_or_likely": sum(
            classifications.get(c, 0)
            for c in ("Benign", "Likely benign")
        ),
        "by_classification": classifications,
        "by_zygosity": zygosities,
        "by_effect": effects,
    }


def _find_igv_snapshots(output_dir: str, variants: List[Dict[str, Any]]) -> Dict[str, str]:
    """IGV 스냅샷 파일을 변이와 매핑합니다."""
    snapshots = {}
    snapshot_dir = os.path.join(output_dir, "snapshots")
    if not os.path.isdir(snapshot_dir):
        # analysis_dir 하위에서도 찾기
        parent = os.path.dirname(output_dir)
        snapshot_dir = os.path.join(parent, "snapshots")

    if not os.path.isdir(snapshot_dir):
        return snapshots

    for img_file in glob.glob(os.path.join(snapshot_dir, "*.png")):
        basename = os.path.basename(img_file)
        # 파일명에서 위치 정보 추출 (예: chr1_12345_A_G.png)
        for v in variants:
            key = f"{v.get('chrom', '')}_{v.get('pos', '')}_{v.get('ref', '')}_{v.get('alt', '')}"
            if key in basename:
                snapshots[v.get("variant_id", "")] = img_file
                break

    return snapshots


# ══════════════════════════════════════════════════════════════
# Variants TSV Export
# ══════════════════════════════════════════════════════════════

TSV_COLUMNS = [
    "variant_id", "chrom", "pos", "ref", "alt",
    "gene", "transcript", "hgvsc", "hgvsp", "effect",
    "zygosity", "gt", "dp", "vaf",
    "gnomad_af", "clinvar_sig_primary", "clinvar_stars",
    "acmg_classification", "acmg_criteria",
    "dbsnp_rsid", "clinvar_dn",
]


def generate_variants_tsv(
    variants: List[Dict[str, Any]],
    output_dir: str,
) -> str:
    """변이 목록을 TSV 파일로 내보냅니다."""
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "variants.tsv")

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for v in variants:
            row = dict(v)
            # acmg_criteria를 문자열로 변환
            criteria = row.get("acmg_criteria", [])
            if isinstance(criteria, list):
                row["acmg_criteria"] = ",".join(criteria)
            writer.writerow(row)

    logger.info(f"Generated variants.tsv: {output_path}")
    return output_path
