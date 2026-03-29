"""
Review Data Generator Module

파이프라인 결과를 기반으로 Portal 리뷰 페이지에 필요한 데이터를 생성합니다.

생성 파일:
    - result.json       : 변이 목록 + annotation + ACMG 분류 (Portal 리뷰 페이지용)
    - qc_summary.json   : QC 메트릭스 요약 (coverage, mapping rate 등)
    - variants.tsv      : 변이 목록 TSV (백업/다운로드용)
    - IGV snapshots     : 주요 변이 위치의 IGV 스냅샷 이미지 (파이프라인에서 생성)

result.json 구조:
    - version, type, generated_at
    - order_id, sample_name, status
    - qc_summary: coverage, alignment, variant_stats
    - variant_stats: total, by_classification, by_zygosity, by_effect
    - variants: 변이 목록 (annotation + ACMG + disease + HPO)
    - disease_panel: BED 기반 질환 패널 정보
    - filter_summary: 적용된 필터 요약
    - igv_snapshots: 변이별 IGV 스냅샷 매핑
    - dark_genes: supplementary unified-pipeline summary (summary/*.txt under analysis_dir);
      optional section_reviews[] (index-aligned with detailed_sections: approved, notes) from Portal;
      customer PDF includes only approved sections
    - metadata: 파이프라인/annotation DB 버전 정보
"""

import os
import re
import csv
import json
import glob
import logging
from typing import Dict, Any, List, Optional, Tuple

from ...datetime_kst import now_kst_iso

logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# QC Summary Extraction
# ══════════════════════════════════════════════════════════════

def _qc_roots(
    analysis_dir: str,
    extra_search_dirs: Optional[List[str]],
    more_search_roots: Optional[List[str]] = None,
) -> List[str]:
    """중복 제거한 QC 검색 루트 (analysis, layout analysis, output, log, multiqc 등)."""
    roots: List[str] = []
    seen = set()

    def add(path: str) -> None:
        if not path or not os.path.isdir(path):
            return
        real = os.path.realpath(path)
        if real in seen:
            return
        seen.add(real)
        roots.append(path)

    add(analysis_dir)
    if extra_search_dirs:
        for d in extra_search_dirs:
            add(d)
    if more_search_roots:
        for d in more_search_roots:
            add(d)
    return roots


def _qc_path_ok(path: str) -> bool:
    """Nextflow work 내 conda env · man 페이지 등 제외."""
    norm = path.replace("\\", "/").lower()
    return "/env/" not in norm and "/viz_env/" not in norm


def _sniff_samtools_flagstat_text(path: str) -> bool:
    """파일명에 flagstat 이 없어도 samtools flagstat 본문이면 True."""
    try:
        if os.path.getsize(path) > 4 * 1024 * 1024:
            return False
    except OSError:
        return False
    try:
        with open(path, "r", errors="ignore") as f:
            chunk = f.read(5000)
    except OSError:
        return False
    low = chunk.lower()
    if "in total" not in low:
        return False
    return (
        "mapped" in low
        or "primary" in low
        or "properly paired" in low
        or re.search(r"\d+\s*\+\s*\d+\s+in total", chunk) is not None
    )


def _collect_extra_flagstat_paths(roots: List[str], already: List[str]) -> List[str]:
    """
    Dark Gene / carrier 레이아웃: multiqc 없이 qc·summary·bam·pipeline_info 아래
    *.txt 등에 flagstat 가 들어 있는 경우.
    """
    seen = set()
    for p in already:
        if os.path.isfile(p):
            try:
                seen.add(os.path.realpath(p))
            except OSError:
                seen.add(p)
    out: List[str] = []
    skip_ext = (
        ".png",
        ".jpg",
        ".jpeg",
        ".gif",
        ".webp",
        ".pdf",
        ".html",
        ".htm",
        ".xml",
        ".bai",
        ".bam",
        ".cram",
        ".fa",
        ".fai",
        ".zip",
        ".bin",
    )
    max_files_per_root = 120
    picked = 0
    subdirs = ("qc", "summary", "bam", "pipeline_info")
    for root in roots:
        picked = 0
        for sd in subdirs:
            base = os.path.join(root, sd)
            if not os.path.isdir(base) or not _qc_path_ok(base):
                continue
            for dirpath, dirnames, filenames in os.walk(base):
                if not _qc_path_ok(dirpath):
                    dirnames[:] = []
                    continue
                rel_sub = os.path.relpath(dirpath, base)
                depth = 0 if rel_sub == "." else rel_sub.count(os.sep) + 1
                if depth >= 12:
                    dirnames[:] = []
                    continue
                for fn in filenames:
                    if picked >= max_files_per_root:
                        break
                    fl = fn.lower()
                    if fl.endswith(skip_ext):
                        continue
                    fp = os.path.join(dirpath, fn)
                    try:
                        if not os.path.isfile(fp):
                            continue
                        rs = os.path.realpath(fp)
                        if rs in seen:
                            continue
                    except OSError:
                        continue
                    if not _sniff_samtools_flagstat_text(fp):
                        continue
                    seen.add(rs)
                    out.append(fp)
                    picked += 1
                if picked >= max_files_per_root:
                    dirnames[:] = []
            if picked >= max_files_per_root:
                break
    return out


def _parse_pipeline_complete_json(path: str) -> Dict[str, Any]:
    """pipeline_complete.json 등에서 정렬·리드 수 키를 휴리스틱으로 수집."""
    result: Dict[str, Any] = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError, TypeError):
        return result

    def to_number(v: Any) -> Optional[float]:
        if isinstance(v, bool):
            return None
        if isinstance(v, (int, float)):
            return float(v)
        if isinstance(v, str):
            t = v.strip().replace(",", "")
            if not t:
                return None
            try:
                return float(t)
            except ValueError:
                return None
        return None

    def consider(k: str, v: Any) -> None:
        fv = to_number(v)
        if fv is None:
            return
        lk = k.lower().replace(" ", "_").replace("-", "_")

        if lk in (
            "total_reads",
            "total_read",
            "reads_total",
            "n_reads",
            "read_pairs",
            "total_fragments",
            "num_reads",
        ):
            result["total_reads"] = int(fv)
        elif lk in (
            "mapped_reads",
            "reads_mapped",
            "n_mapped",
            "aligned_reads",
            "reads_aligned",
        ):
            result["mapped_reads"] = int(fv)
        elif lk in (
            "duplicates",
            "duplicate_reads",
            "read_duplicates",
            "duplicate_read_pairs",
        ):
            result["duplicates"] = int(fv)
        elif lk in (
            "properly_paired",
            "reads_properly_paired",
            "proper_pairs",
            "properly_paired_reads",
        ):
            result["properly_paired"] = int(fv)
        elif (
            "mapping" in lk or lk.startswith("pct_mapped") or "percent_mapped" in lk
        ) and ("unmapped" not in lk):
            if 0 < fv <= 1.0:
                result["mapping_rate"] = round(100.0 * fv, 4)
            elif fv > 1.0:
                result["mapping_rate"] = round(fv, 4)
        elif "properly" in lk and "paired" in lk and (
            "pct" in lk or "rate" in lk or "percent" in lk or lk.endswith("_pct")
        ):
            if 0 < fv <= 1.0:
                result["properly_paired_rate"] = round(100.0 * fv, 4)
            elif fv > 1.0:
                result["properly_paired_rate"] = round(fv, 4)

    def walk(obj: Any) -> None:
        if isinstance(obj, dict):
            for k, v in obj.items():
                consider(k, v)
                walk(v)
        elif isinstance(obj, list):
            for it in obj[:400]:
                walk(it)

    walk(data)

    tr = result.get("total_reads")
    mr = result.get("mapped_reads")
    if tr and mr is not None and result.get("mapping_rate") is None and tr > 0:
        result["mapping_rate"] = round(100.0 * float(mr) / float(tr), 3)
    if (
        tr
        and result.get("properly_paired") is not None
        and result.get("properly_paired_rate") is None
        and tr > 0
    ):
        result["properly_paired_rate"] = round(
            100.0 * float(result["properly_paired"]) / float(tr), 3
        )
    return result


def _split_metric_kv_line(line: str) -> Optional[tuple]:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    if "\t" in line:
        parts = line.split("\t", 1)
        if len(parts) == 2 and parts[0].strip():
            return parts[0].strip(), parts[1].strip()
    if "=" in line and not line.startswith("="):
        k, _, v = line.partition("=")
        if k.strip() and v.strip():
            return k.strip(), v.strip()
    if ":" in line:
        k, _, v = line.partition(":")
        if (
            k.strip()
            and v.strip()
            and len(k) < 96
            and "://" not in line
            and not re.match(r"^\d{4}-\d{2}-\d{2}", k.strip())
        ):
            return k.strip(), v.strip()
    return None


def _apply_flat_qc_metric(
    key: str,
    val: str,
    alignment: Dict[str, Any],
    coverage: Dict[str, Any],
) -> None:
    """qc/*_qc_metrics.txt 등 단일 키–값(정규화된 key)."""
    key = key.lower().replace(" ", "_").replace("-", "_")
    raw = val.strip().strip('"').strip("'")
    pct_suffix = raw.endswith("%")
    num_s = raw.rstrip("%").replace(",", "").strip()
    try:
        fv = float(num_s)
    except ValueError:
        return

    def as_int(x: float) -> int:
        return int(round(x))

    if any(
        x in key
        for x in (
            "mean_coverage",
            "avg_coverage",
            "average_coverage",
            "mean_depth",
            "avg_depth",
            "target_mean_coverage",
        )
    ):
        coverage["mean_coverage"] = fv
    elif "min_coverage" in key or key.endswith("coverage_min"):
        coverage["min_coverage"] = fv
    elif "max_coverage" in key or key.endswith("coverage_max"):
        coverage["max_coverage"] = fv
    elif "on_target" in key or "ontarget" in key:
        if fv <= 1.0:
            coverage["on_target_rate"] = round(100.0 * fv, 4)
        else:
            coverage["on_target_rate"] = fv
    elif (
        any(
            x in key
            for x in (
                "total_read",
                "total_sequences",
                "total_fragments",
                "read_pairs",
                "read_pair_count",
            )
        )
        and "rate" not in key
        and "percent" not in key
        and "pct" not in key
        and "mapped" not in key
    ):
        alignment["total_reads"] = as_int(fv)
    elif (
        key in ("mapped_reads", "reads_mapped", "aligned_reads", "n_mapped")
        or ("mapped" in key and "read" in key and "unmapped" not in key)
    ) and "rate" not in key and "percent" not in key and "pct" not in key and not pct_suffix:
        alignment["mapped_reads"] = as_int(fv)
    elif "duplicate" in key and "read" in key and "rate" not in key:
        alignment["duplicates"] = as_int(fv)
    elif "properly_paired" in key and "rate" not in key and "percent" not in key and "pct" not in key:
        alignment["properly_paired"] = as_int(fv)
    elif (
        "mapping" in key and ("rate" in key or "percent" in key or "pct" in key)
    ) or key in ("map_rate", "mapped_pct", "percent_mapped", "pct_mapped"):
        if fv <= 1.0 and not pct_suffix:
            alignment["mapping_rate"] = round(100.0 * fv, 4)
        else:
            alignment["mapping_rate"] = round(min(fv, 100.0), 4)
    elif (
        "properly" in key
        and "paired" in key
        and ("rate" in key or "percent" in key or "pct" in key or pct_suffix)
    ):
        if fv <= 1.0 and not pct_suffix:
            alignment["properly_paired_rate"] = round(100.0 * fv, 4)
        else:
            alignment["properly_paired_rate"] = round(min(fv, 100.0), 4)
    elif ("average" in key and "quality" in key) or key in ("mean_quality", "avg_q", "mean_q"):
        alignment["average_quality"] = fv
    else:
        # % Bases > Nx  /  ≥ Nx  (e.g. "% Bases > 20x", "Pct_Bases_>=_50x")
        m_depth = re.search(r"(\d+)x", key)
        if m_depth and ("base" in key or "cov" in key or "depth" in key):
            threshold = m_depth.group(1)
            if pct_suffix:
                coverage[f"pct_bases_{threshold}x"] = fv
            elif fv <= 1.0:
                coverage[f"pct_bases_{threshold}x"] = round(100.0 * fv, 2)
            else:
                coverage[f"pct_bases_{threshold}x"] = fv


def _parse_carrier_qc_flat_txt(path: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    qc/*_qc_metrics.txt, *target_coverage*.txt
    TSV 헤더+1행 또는 key:value / key=value / 탭 구분 행.
    Returns (alignment_dict, coverage_dict).
    """
    alignment: Dict[str, Any] = {}
    coverage: Dict[str, Any] = {}
    try:
        with open(path, "r", errors="replace") as f:
            raw_lines = f.readlines()
    except OSError:
        return alignment, coverage
    lines = [ln.strip() for ln in raw_lines if ln.strip() and not ln.strip().startswith("#")]
    if len(lines) >= 2 and "\t" in lines[0]:
        hdr = [x.strip() for x in lines[0].split("\t")]
        if len(hdr) >= 2 and "\t" in lines[1]:
            vals = [x.strip() for x in lines[1].split("\t")]
            if len(hdr) == len(vals):
                all_headers = all(
                    re.match(r"^[A-Za-z_][A-Za-z0-9_ %()]*$", h) for h in hdr
                )
                if all_headers:
                    for k, v in zip(hdr, vals):
                        kn = k.lower().strip().replace(" ", "_").replace("-", "_")
                        kn = re.sub(r"[%()]+", "", kn)
                        _apply_flat_qc_metric(kn, v, alignment, coverage)
                    return alignment, coverage
    for ln in lines:
        pair = _split_metric_kv_line(ln)
        if not pair:
            continue
        k, v = pair
        kn = k.lower().strip().replace(" ", "_").replace("-", "_")
        kn = re.sub(r"[%()]+", "", kn)
        _apply_flat_qc_metric(kn, v, alignment, coverage)
    return alignment, coverage


def _merge_alignment_dicts(dicts: List[Dict[str, Any]]) -> Dict[str, Any]:
    """여러 flagstat/stats 파싱 결과를 합칩니다. total_reads가 큰 결과를 베이스로 하고 나머지로 빈 키만 채웁니다."""
    non_empty = [d for d in dicts if d]
    if not non_empty:
        return {}
    scored = sorted(
        non_empty,
        key=lambda d: (
            int(d.get("total_reads") or d.get("raw_total_sequences") or 0),
            sum(1 for v in d.values() if v is not None and v != ""),
        ),
        reverse=True,
    )
    out: Dict[str, Any] = dict(scored[0])
    for d in scored[1:]:
        for k, v in d.items():
            if k not in out or out[k] is None:
                out[k] = v
    return out


def _alignment_postprocess(aln: Dict[str, Any]) -> None:
    """alignment 맵을 보강 (파생 필드·이름 정규화)."""
    if not aln:
        return
    if aln.get("total_reads") in (None, 0) and aln.get("raw_total_sequences"):
        aln["total_reads"] = aln["raw_total_sequences"]
    rt = aln.get("total_reads") or aln.get("raw_total_sequences")
    mr = aln.get("mapped_reads")
    if rt and mr is not None and aln.get("mapping_rate") is None:
        try:
            rt_f = float(rt)
            if rt_f > 0:
                aln["mapping_rate"] = round(100.0 * float(mr) / rt_f, 3)
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    pp = aln.get("properly_paired")
    if (
        rt
        and pp is not None
        and aln.get("properly_paired_rate") is None
    ):
        try:
            rt_f = float(rt)
            if rt_f > 0:
                aln["properly_paired_rate"] = round(100.0 * float(pp) / rt_f, 3)
        except (TypeError, ValueError, ZeroDivisionError):
            pass


def _peek_is_samtools_stats(path: str) -> bool:
    try:
        with open(path, "r", errors="ignore") as f:
            chunk = f.read(4096)
    except OSError:
        return False
    return "This file was produced by samtools stats" in chunk or "\nSN\t" in chunk


def extract_qc_summary(
    analysis_dir: str,
    sample_name: str,
    extra_search_dirs: Optional[List[str]] = None,
    more_search_roots: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    파이프라인 분석 디렉토리에서 QC 메트릭스를 추출합니다.

    extra_search_dirs: artifact 전용 경로(서비스 데몬 work 루트)와 달리,
    Nextflow가 mosdepth 등을 쓰는 `layout_base/analysis/{work}/{sample}/` 를 추가로 지정합니다.
    more_search_roots: output_dir, layout_base/output, log 등 (MultiQC·Picard·publish 디렉터리).

    주요 QC 파일:
        - *mosdepth*summary*  (coverage)
        - *flagstat* (파일명에 flagstat 포함; .flagstat.txt 등)
        - samtools stats: *{sample_name}*.stats, *.bam.stats, 또는 헤더가 samtools.stats
        - multiqc_data/multiqc_general_stats.txt, mqc_samtools_flagstat*.txt
        - Picard *alignment_summary_metrics*
        - *coverage_report*
        - qc/*_qc_metrics.txt, qc/*target_coverage*.txt (Dark Gene / Twist Exome 스타일 플랫 QC)
    """
    qc = {
        "sample_name": sample_name,
        "generated_at": now_kst_iso(),
        "coverage": {},
        "alignment": {},
        "variant_stats": {},
    }

    roots = _qc_roots(analysis_dir, extra_search_dirs, more_search_roots)

    for root in roots:
        for mf in glob.glob(os.path.join(root, "**", "*mosdepth*summary*"), recursive=True):
            if not _qc_path_ok(mf):
                continue
            try:
                qc["coverage"].update(_parse_mosdepth_summary(mf))
            except Exception as e:
                logger.warning(f"Failed to parse mosdepth file {mf}: {e}")

    flagstat_paths: List[str] = []
    for root in roots:
        for ff in glob.glob(os.path.join(root, "**", "*flagstat*"), recursive=True):
            if not os.path.isfile(ff) or not _qc_path_ok(ff):
                continue
            low = ff.lower()
            if low.endswith((".md", ".html", ".xml", ".json")):
                continue
            flagstat_paths.append(ff)

    flagstat_paths.extend(_collect_extra_flagstat_paths(roots, flagstat_paths))

    flagstat_parsed: List[Dict[str, Any]] = []
    for ff in flagstat_paths:
        try:
            p = _parse_flagstat(ff)
            if p:
                flagstat_parsed.append(p)
        except Exception as e:
            logger.warning(f"Failed to parse flagstat file {ff}: {e}")

    stats_paths: List[str] = []
    stats_globs = [
        f"*{sample_name}*.stats",
        f"*{sample_name}*.stats.txt",
        "*samtools*stats*.txt",
        "*.bam.stats",
    ]
    for root in roots:
        for pattern in stats_globs:
            for sf in glob.glob(os.path.join(root, "**", pattern), recursive=True):
                if not os.path.isfile(sf) or not _qc_path_ok(sf):
                    continue
                if sf not in stats_paths:
                    stats_paths.append(sf)
        extra_stats: List[str] = []
        for pattern in ("*.stats", "*.stats.txt"):
            for sf in glob.glob(os.path.join(root, "**", pattern), recursive=True):
                if not os.path.isfile(sf) or not _qc_path_ok(sf) or sf in stats_paths:
                    continue
                try:
                    if os.path.getsize(sf) > 25 * 1024 * 1024:
                        continue
                except OSError:
                    continue
                extra_stats.append(sf)
        norm_s = sample_name.replace("\\", "/")
        extra_stats.sort(
            key=lambda p: (
                0 if sample_name in os.path.basename(p) else 1,
                0 if norm_s in p.replace("\\", "/") else 1,
                len(p),
            )
        )
        for sf in extra_stats[:120]:
            if _peek_is_samtools_stats(sf):
                stats_paths.append(sf)

    stats_parsed: List[Dict[str, Any]] = []
    for sf in stats_paths:
        try:
            parsed = _parse_samtools_stats(sf)
            if parsed:
                stats_parsed.append(parsed)
        except Exception as e:
            logger.warning(f"Failed to parse stats file {sf}: {e}")

    star_parsed: List[Dict[str, Any]] = []
    qualimap_parsed: List[Dict[str, Any]] = []
    fastp_parsed: List[Dict[str, Any]] = []
    multiqc_parsed: List[Dict[str, Any]] = []
    picard_parsed: List[Dict[str, Any]] = []
    pipeline_complete_parsed: List[Dict[str, Any]] = []
    for root in roots:
        for lf in glob.glob(os.path.join(root, "**", "Log.final.out"), recursive=True):
            if not os.path.isfile(lf) or not _qc_path_ok(lf):
                continue
            try:
                if os.path.getsize(lf) > 8 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                if not _peek_is_star_log_final(lf):
                    continue
                p = _parse_star_log_final(lf)
                if p:
                    star_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse STAR Log.final.out {lf}: {e}")
        for gf in glob.glob(
            os.path.join(root, "**", "genome_results.txt"), recursive=True
        ):
            if not os.path.isfile(gf) or not _qc_path_ok(gf):
                continue
            try:
                if os.path.getsize(gf) > 2 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                if not _peek_is_qualimap_genome_results(gf):
                    continue
                p = _parse_qualimap_genome_results(gf)
                if p:
                    qualimap_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse QualiMap genome_results {gf}: {e}")
        for jp in glob.glob(os.path.join(root, "**", "*fastp*.json"), recursive=True):
            if not os.path.isfile(jp) or not _qc_path_ok(jp):
                continue
            low = jp.lower()
            if "multiqc" in low and "data" in low:
                continue
            try:
                if os.path.getsize(jp) > 12 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                p = _parse_fastp_json(jp)
                if p:
                    fastp_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse fastp json {jp}: {e}")
        for mq in glob.glob(
            os.path.join(root, "**", "multiqc_general_stats.txt"), recursive=True
        ):
            if not os.path.isfile(mq) or not _qc_path_ok(mq):
                continue
            try:
                if os.path.getsize(mq) > 6 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                p = _parse_multiqc_tsv_alignment(mq, sample_name)
                if p:
                    multiqc_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse MultiQC general stats {mq}: {e}")
        for mq in glob.glob(
            os.path.join(root, "**", "mqc_samtools_flagstat*.txt"), recursive=True
        ):
            if not os.path.isfile(mq) or not _qc_path_ok(mq):
                continue
            low = mq.lower()
            if low.endswith(".md"):
                continue
            try:
                if os.path.getsize(mq) > 2 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                p = _parse_multiqc_tsv_alignment(mq, sample_name)
                if p:
                    multiqc_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse MultiQC samtools flagstat {mq}: {e}")
        for pp in glob.glob(
            os.path.join(root, "**", "*alignment_summary_metrics*"), recursive=True
        ):
            if not os.path.isfile(pp) or not _qc_path_ok(pp):
                continue
            low = pp.lower()
            if low.endswith((".md", ".html", ".pdf")):
                continue
            try:
                if os.path.getsize(pp) > 5 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                if not _peek_picard_alignment_summary(pp):
                    continue
                p = _parse_picard_alignment_summary_metrics(pp)
                if p:
                    picard_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse Picard alignment metrics {pp}: {e}")

        # Picard MarkDuplicates  →  total_reads, mapping_rate, duplicates
        for dp in glob.glob(
            os.path.join(root, "**", "*duplicate_metrics*"), recursive=True
        ):
            if not os.path.isfile(dp) or not _qc_path_ok(dp):
                continue
            if dp.lower().endswith((".md", ".html", ".pdf")):
                continue
            try:
                if os.path.getsize(dp) > 5 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                p = _parse_picard_duplicate_metrics(dp)
                if p:
                    picard_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse Picard duplicate metrics {dp}: {e}")

        # samtools idxstats  →  total_reads, mapped_reads, mapping_rate
        for ix in glob.glob(
            os.path.join(root, "**", "*idxstats*"), recursive=True
        ):
            if not os.path.isfile(ix) or not _qc_path_ok(ix):
                continue
            if ix.lower().endswith((".md", ".html", ".pdf", ".bai", ".bam", ".cram")):
                continue
            try:
                if os.path.getsize(ix) > 2 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                p = _parse_samtools_idxstats(ix)
                if p:
                    picard_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse idxstats {ix}: {e}")
        for pc in glob.glob(
            os.path.join(root, "**", "pipeline_complete.json"), recursive=True
        ):
            if not os.path.isfile(pc) or not _qc_path_ok(pc):
                continue
            try:
                if os.path.getsize(pc) > 20 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                p = _parse_pipeline_complete_json(pc)
                if p:
                    pipeline_complete_parsed.append(p)
            except Exception as e:
                logger.warning(f"Failed to parse pipeline_complete.json {pc}: {e}")

    carrier_flat_parsed: List[Dict[str, Any]] = []
    for root in roots:
        for qf in glob.glob(os.path.join(root, "**", "qc", "*.txt"), recursive=True):
            if not os.path.isfile(qf) or not _qc_path_ok(qf):
                continue
            bn = os.path.basename(qf).lower()
            if "qc_metrics" not in bn and "target_coverage" not in bn:
                continue
            try:
                if os.path.getsize(qf) > 3 * 1024 * 1024:
                    continue
            except OSError:
                continue
            try:
                aln, cov = _parse_carrier_qc_flat_txt(qf)
                if cov:
                    qc["coverage"].update(cov)
                if aln:
                    carrier_flat_parsed.append(aln)
            except Exception as e:
                logger.warning(f"Failed to parse carrier qc flat txt {qf}: {e}")

    alignment_parts: List[Dict[str, Any]] = (
        flagstat_parsed
        + stats_parsed
        + star_parsed
        + qualimap_parsed
        + fastp_parsed
        + multiqc_parsed
        + picard_parsed
        + pipeline_complete_parsed
        + carrier_flat_parsed
    )
    qc["alignment"] = _merge_alignment_dicts(alignment_parts)
    _alignment_postprocess(qc["alignment"])

    for root in roots:
        for cf in glob.glob(os.path.join(root, "**", "*coverage_report*"), recursive=True):
            if not _qc_path_ok(cf):
                continue
            try:
                qc["coverage"].update(_parse_coverage_report(cf))
            except Exception as e:
                logger.warning(f"Failed to parse coverage report {cf}: {e}")

    if qc["coverage"] and not qc["alignment"]:
        logger.info(
            "extract_qc_summary: no alignment metrics parsed for sample=%s (%d search roots); "
            "look for qc/* text (flagstat), pipeline_complete.json, or Picard under output",
            sample_name,
            len(roots),
        )

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
    """samtools / sambamba flagstat 파싱 (형식·확장자 변형 허용)."""
    result: Dict[str, Any] = {}
    with open(path, "r", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            low = line.lower()
            if "in total" in low:
                m = re.search(r"(\d+)\s*\+\s*\d+\s+in total\b", line)
                if not m:
                    # 일부 로그: "12345 in total (QC-passed reads + QC-failed reads)"
                    m = re.search(r"^(\d+)\s+in total\b", line)
                if m and result.get("total_reads") is None:
                    result["total_reads"] = int(m.group(1))
                continue
            if "mapped" in low and "%" in line and "unmapped" not in low:
                m = re.search(
                    r"(\d+)\s*\+\s*\d+\s+mapped\s+\(([\d.]+)\s*%",
                    line,
                    re.IGNORECASE,
                )
                if m:
                    result["mapped_reads"] = int(m.group(1))
                    result["mapping_rate"] = float(m.group(2))
                continue
            if "properly paired" in low and "%" in line:
                m = re.search(
                    r"(\d+)\s*\+\s*\d+\s+properly paired\s+\(([\d.]+)\s*%",
                    line,
                    re.IGNORECASE,
                )
                if m:
                    result["properly_paired"] = int(m.group(1))
                    result["properly_paired_rate"] = float(m.group(2))
                continue
            if "duplicate" in low:
                m = re.search(r"(\d+)\s*\+\s*\d+\s+duplicates?\b", line, re.IGNORECASE)
                if m:
                    result["duplicates"] = int(m.group(1))
    return result


def _peek_is_star_log_final(path: str) -> bool:
    try:
        with open(path, "r", errors="ignore") as f:
            head = f.read(2500)
    except OSError:
        return False
    return "STAR" in head and "Number of input reads" in head


def _parse_star_log_final(path: str) -> Dict[str, Any]:
    """STAR aligner Log.final.out 에서 총 리드·매핑( unique + multi ) 추출."""
    result: Dict[str, Any] = {}
    try:
        with open(path, "r", errors="replace") as f:
            text = f.read()
    except OSError:
        return result

    def grab_int(label_pat: str) -> Optional[int]:
        m = re.search(label_pat, text, re.MULTILINE)
        if not m:
            return None
        try:
            return int(m.group(1))
        except (ValueError, IndexError):
            return None

    tr = grab_int(r"Number of input reads\s*\|\s*(\d+)")
    if tr is not None:
        result["total_reads"] = tr
    um = grab_int(r"Uniquely mapped reads number\s*\|\s*(\d+)")
    mm = grab_int(r"Number of reads mapped to multiple loci\s*\|\s*(\d+)")
    if um is not None:
        mm_i = mm if mm is not None else 0
        result["mapped_reads"] = um + mm_i
        if tr:
            result["mapping_rate"] = round(100.0 * (um + mm_i) / tr, 3)
    return result


def _peek_is_qualimap_genome_results(path: str) -> bool:
    try:
        with open(path, "r", errors="ignore") as f:
            head = f.read(1200)
    except OSError:
        return False
    return "QualiMap" in head or (
        ">>>>>>>" in head and "Number of reads" in head
    )


def _parse_qualimap_genome_results(path: str) -> Dict[str, Any]:
    """QualiMap genome_results.txt (reads / mapped / % 일부 변형)."""
    result: Dict[str, Any] = {}
    try:
        with open(path, "r", errors="replace") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                low = line.lower()
                m_nr = re.match(
                    r"number of reads\s*[:=]\s*([\d,]+)", line, re.IGNORECASE
                )
                if m_nr:
                    result["total_reads"] = int(m_nr.group(1).replace(",", ""))
                    continue
                m_rm = re.match(r"reads mapped\s*:\s*([\d,]+)", line, re.IGNORECASE)
                if m_rm:
                    result["mapped_reads"] = int(m_rm.group(1).replace(",", ""))
                    continue
                if "mapped" in low and "reads" in low and "%" in line:
                    m_pct = re.search(r"([\d.]+)\s*%", line)
                    if m_pct and result.get("mapping_rate") is None:
                        result["mapping_rate"] = float(m_pct.group(1))
    except OSError:
        return result
    tr = result.get("total_reads")
    mr = result.get("mapped_reads")
    if tr and mr is not None and result.get("mapping_rate") is None and tr > 0:
        result["mapping_rate"] = round(100.0 * float(mr) / float(tr), 3)
    return result


def _parse_fastp_json(path: str) -> Dict[str, Any]:
    """fastp JSON summary 에서 총 리드·중복률·평균 Q(가능 시)."""
    result: Dict[str, Any] = {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError, TypeError):
        return result
    if not isinstance(data, dict):
        return result
    summ = data.get("summary") or {}
    bf = summ.get("before_filtering") or {}
    tr = bf.get("total_reads")
    if tr is not None:
        try:
            result["total_reads"] = int(tr)
        except (TypeError, ValueError):
            pass
    dup = data.get("duplication") or {}
    r = dup.get("rate")
    if r is not None and result.get("total_reads"):
        try:
            rv = float(r)
            if 0 <= rv <= 1:
                result["duplicates"] = int(round(result["total_reads"] * rv))
            elif rv > 1:
                result["duplicates"] = int(round(rv))
        except (TypeError, ValueError):
            pass
    # before_filtering average read quality if present
    for key in ("read1_mean_quality", "quality_mean", "mean_quality"):
        q = bf.get(key)
        if q is not None:
            try:
                result["average_quality"] = float(q)
                break
            except (TypeError, ValueError):
                pass
    return result


def _parse_int_cell(s: str) -> Optional[int]:
    t = (s or "").strip().replace(",", "").replace("%", "")
    if not t or t.lower() in ("na", "n/a", "."):
        return None
    try:
        return int(float(t))
    except ValueError:
        return None


def _parse_float_cell(s: str) -> Optional[float]:
    t = (s or "").strip().replace(",", "").replace("%", "")
    if not t or t.lower() in ("na", "n/a", "."):
        return None
    try:
        return float(t)
    except ValueError:
        return None


def _parse_multiqc_tsv_alignment(path: str, sample_name: str) -> Dict[str, Any]:
    """
    MultiQC 탭 표 (multiqc_general_stats.txt, mqc_samtools_flagstat*.txt 등).
    헤더 이름에서 total/mapped/properly paired/duplicate 열을 추정합니다.
    """
    result: Dict[str, Any] = {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            rows_raw: List[List[str]] = []
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rows_raw.append(line.split("\t"))
    except OSError:
        return result
    if len(rows_raw) < 2:
        return result
    header = [c.strip() for c in rows_raw[0]]
    body = rows_raw[1:]
    sample_col = 0 if header and header[0].lower() == "sample" else None

    ti = mi = mri = ppi = ppri = dupi = None
    for i, h in enumerate(header):
        hl = h.lower().replace(" ", "_")
        if sample_col is not None and i == sample_col:
            continue
        if "unmapped" in hl:
            continue
        if "total" in hl and "read" in hl and "base" not in hl and ti is None:
            ti = i
        if "mapped" in hl and "read" in hl:
            if any(x in hl for x in ("percent", "pct", "_ratio")) or "%" in h:
                if mri is None:
                    mri = i
            elif mi is None:
                mi = i
        if "properly" in hl and "paired" in hl:
            if any(x in hl for x in ("percent", "pct")) or "%" in h:
                if ppri is None:
                    ppri = i
            elif ppi is None:
                ppi = i
        if "duplicate" in hl and any(x in hl for x in ("read", "fragment")):
            if dupi is None and "percent" not in hl:
                dupi = i

    candidates = body
    if sample_col is not None and sample_name and body:
        sn = sample_name.strip()
        matched = [
            r
            for r in body
            if len(r) > sample_col
            and (
                r[sample_col].strip() == sn
                or sn in r[sample_col]
                or r[sample_col].strip() in sn
            )
        ]
        if matched:
            candidates = matched

    chosen: Optional[List[str]] = None
    if ti is not None and candidates:
        best_tot = -1
        for row in candidates:
            if len(row) <= ti:
                continue
            t = _parse_int_cell(row[ti])
            if t is not None and t > best_tot:
                best_tot = t
                chosen = row
    if chosen is None and candidates:
        chosen = candidates[0]
    if chosen is None or len(chosen) < len(header):
        return result

    if ti is not None:
        v = _parse_int_cell(chosen[ti])
        if v is not None:
            result["total_reads"] = v
    if mi is not None:
        v = _parse_int_cell(chosen[mi])
        if v is not None:
            result["mapped_reads"] = v
    if mri is not None:
        v = _parse_float_cell(chosen[mri])
        if v is not None:
            result["mapping_rate"] = v
    if ppi is not None:
        v = _parse_int_cell(chosen[ppi])
        if v is not None:
            result["properly_paired"] = v
    if ppri is not None:
        v = _parse_float_cell(chosen[ppri])
        if v is not None:
            result["properly_paired_rate"] = v
    if dupi is not None:
        v = _parse_int_cell(chosen[dupi])
        if v is not None:
            result["duplicates"] = v

    tr = result.get("total_reads")
    mr = result.get("mapped_reads")
    if tr and mr is not None and result.get("mapping_rate") is None and tr > 0:
        result["mapping_rate"] = round(100.0 * float(mr) / float(tr), 3)
    if (
        tr
        and result.get("properly_paired") is not None
        and result.get("properly_paired_rate") is None
        and tr > 0
    ):
        result["properly_paired_rate"] = round(
            100.0 * float(result["properly_paired"]) / float(tr), 3
        )
    return result


def _parse_picard_duplicate_metrics(path: str) -> Dict[str, Any]:
    """
    Picard MarkDuplicates *duplicate_metrics* 파일 파싱.
    LIBRARY 행에서 total_reads, mapped_reads, mapping_rate, duplicates, duplicate_rate 추출.
    """
    result: Dict[str, Any] = {}
    try:
        with open(path, "r", errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return result

    header: List[str] = []
    for i, line in enumerate(lines):
        if line.startswith("LIBRARY\t"):
            header = line.strip().split("\t")
            # first data row
            for data_line in lines[i + 1:]:
                data_line = data_line.strip()
                if not data_line or data_line.startswith("#") or data_line.startswith("BIN"):
                    break
                parts = data_line.split("\t")
                if len(parts) < len(header):
                    break
                idx = {name: j for j, name in enumerate(header)}

                def _int(col: str) -> Optional[int]:
                    j = idx.get(col)
                    return _parse_int_cell(parts[j]) if j is not None else None

                def _flt(col: str) -> Optional[float]:
                    j = idx.get(col)
                    return _parse_float_cell(parts[j]) if j is not None else None

                pairs = _int("READ_PAIRS_EXAMINED") or 0
                unpaired = _int("UNPAIRED_READS_EXAMINED") or 0
                unmapped = _int("UNMAPPED_READS") or 0
                secondary = _int("SECONDARY_OR_SUPPLEMENTARY_RDS") or 0
                pair_dups = _int("READ_PAIR_DUPLICATES") or 0
                unpaired_dups = _int("UNPAIRED_READ_DUPLICATES") or 0
                pct_dup = _flt("PERCENT_DUPLICATION")

                # Primary reads only (exclude secondary/supp)
                total = pairs * 2 + unpaired + unmapped
                if total > 0:
                    result["total_reads"] = total
                    result["mapped_reads"] = total - unmapped
                    result["mapping_rate"] = round((total - unmapped) / total * 100, 4)

                dups = pair_dups * 2 + unpaired_dups
                if dups > 0:
                    result["duplicates"] = dups
                if pct_dup is not None:
                    result["duplicate_rate"] = round(pct_dup * 100, 4)
                break
            break
    return result


def _parse_samtools_idxstats(path: str) -> Dict[str, Any]:
    """
    samtools idxstats: chrom\\tlength\\tmapped\\tunmapped 형식.
    * (unplaced/unmapped) 행까지 합산하여 total_reads / mapping_rate 계산.
    """
    result: Dict[str, Any] = {}
    total_mapped = 0
    total_unmapped = 0
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                try:
                    total_mapped += int(parts[2])
                    total_unmapped += int(parts[3])
                except ValueError:
                    continue
    except OSError:
        return result

    total = total_mapped + total_unmapped
    if total > 0:
        result["total_reads"] = total
        result["mapped_reads"] = total_mapped
        result["mapping_rate"] = round(total_mapped / total * 100, 4)
    return result


def _peek_picard_alignment_summary(path: str) -> bool:
    try:
        with open(path, "r", errors="ignore") as f:
            head = f.read(6000)
    except OSError:
        return False
    return (
        "AlignmentSummaryMetrics" in head
        and "CATEGORY" in head
        and "TOTAL_READS" in head
    )


def _parse_picard_alignment_summary_metrics(path: str) -> Dict[str, Any]:
    """Picard CollectAlignmentSummaryMetrics (PAIR / UNPAIRED 행)."""
    result: Dict[str, Any] = {}
    try:
        with open(path, "r", errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return result
    header: List[str] = []
    start = -1
    for i, line in enumerate(lines):
        if line.startswith("CATEGORY\t"):
            header = line.strip().split("\t")
            start = i
            break
    if start < 0 or not header:
        return result
    idx = {name: j for j, name in enumerate(header)}

    def pull_row(parts: List[str]) -> None:
        if len(parts) < len(header):
            return
        tr = idx.get("TOTAL_READS")
        if tr is not None:
            v = _parse_int_cell(parts[tr])
            if v is not None:
                result["total_reads"] = v
        am = idx.get("PF_READS_ALIGNED")
        if am is not None:
            v = _parse_int_cell(parts[am])
            if v is not None:
                result["mapped_reads"] = v
        pct = idx.get("PF_READS_ALIGNED_PCT")
        if pct is not None:
            v = _parse_float_cell(parts[pct])
            if v is not None:
                result["mapping_rate"] = v
        pp = idx.get("PF_READS_PROPERLY_PAIRED")
        if pp is not None:
            v = _parse_int_cell(parts[pp])
            if v is not None:
                result["properly_paired"] = v
        ppc = idx.get("PF_READS_PROPERLY_PAIRED_PCT")
        if ppc is not None:
            v = _parse_float_cell(parts[ppc])
            if v is not None:
                result["properly_paired_rate"] = v

    for line in lines[start + 1 :]:
        if line.startswith("PAIR\t") or line.startswith("UNPAIRED\t"):
            pull_row(line.strip().split("\t"))
            if result.get("total_reads"):
                break
    return result


def _parse_samtools_stats(path: str) -> Dict[str, Any]:
    """samtools stats 파싱 (SN 테이블)."""
    result: Dict[str, Any] = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            if not line.startswith("SN\t"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            key = parts[1].rstrip(":")
            val = parts[2]
            try:
                if key == "raw total sequences":
                    result["raw_total_sequences"] = int(val)
                elif key == "sequences":
                    result.setdefault("raw_total_sequences", int(val))
                elif key == "reads mapped":
                    result.setdefault("mapped_reads", int(val))
                elif key == "reads duplicated":
                    result["duplicates"] = int(val)
                elif key == "reads properly paired":
                    result["properly_paired"] = int(val)
                elif key in (
                    "reads properly paired %",
                    "reads properly paired (%)",
                    "percentage of properly paired reads",
                ):
                    result["properly_paired_rate"] = float(val)
                elif key == "insert size average":
                    result["insert_size_avg"] = float(val)
                elif key == "insert size standard deviation":
                    result["insert_size_std"] = float(val)
                elif key == "average quality":
                    result["average_quality"] = float(val)
            except (ValueError, TypeError):
                continue
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
# QC metric chart images (qc/, summary/ under output roots)
# ══════════════════════════════════════════════════════════════

_QC_IMAGE_EXT = (".png", ".jpg", ".jpeg", ".svg", ".webp", ".gif")


def _metric_keys_for_qc_image_filename(basename_lower: str) -> List[str]:
    """파일명 토큰으로 qc_summary coverage/alignment 키에 연결."""
    keys: List[str] = []
    bn = basename_lower

    if any(x in bn for x in ("qual", "quality", "qscore")):
        keys.extend(["average_quality"])
    if "insert" in bn:
        keys.extend(["insert_size_avg", "insert_size_std"])
    if "dup" in bn or "duplicate" in bn:
        keys.extend(["duplicates", "duplicate_rate"])
    if "proper" in bn or ("paired" in bn and "unpair" not in bn):
        keys.extend(["properly_paired", "properly_paired_rate"])
    if any(x in bn for x in ("map", "align", "flagstat")):
        keys.extend(["mapping_rate", "mapped_reads"])
    if ("read" in bn or "total" in bn) and "qual" not in bn:
        if any(x in bn for x in ("total", "count", "seq")):
            keys.extend(["total_reads", "raw_total_sequences"])
    if any(
        x in bn
        for x in (
            "mosdepth",
            "depth",
            "coverage",
            "target",
            "on_target",
            "fraction",
            "breadth",
        )
    ):
        keys.extend(
            ["mean_coverage", "min_coverage", "max_coverage", "on_target_rate"]
        )
    for m in re.finditer(r"(\d+)x", bn):
        keys.append(f"pct_bases_{m.group(1)}x")

    out: List[str] = []
    seen = set()
    for k in keys:
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


def discover_qc_metric_images(search_roots: List[str]) -> Dict[str, List[str]]:
    """
    output/layout 루트 아래 qc/, summary/ 에서 이미지를 찾아 메트릭 키별로 그룹화.
    값은 job.output_dir 기준 상대 경로 (예: qc/depth.png) — download API filename 으로 사용.
    """
    metric_to_paths: Dict[str, List[str]] = {}
    roots_done: set = set()
    for root in search_roots:
        if not root or not os.path.isdir(root):
            continue
        try:
            rp = os.path.realpath(root)
        except OSError:
            continue
        if rp in roots_done:
            continue
        roots_done.add(rp)
        for sub in ("qc", "summary"):
            d = os.path.join(root, sub)
            if not os.path.isdir(d):
                continue
            try:
                for ent in sorted(os.listdir(d)):
                    if ent.startswith("."):
                        continue
                    path = os.path.join(d, ent)
                    if not os.path.isfile(path):
                        continue
                    low = ent.lower()
                    if not low.endswith(_QC_IMAGE_EXT):
                        continue
                    rel = f"{sub}/{ent}".replace("\\", "/")
                    mkeys = _metric_keys_for_qc_image_filename(low)
                    for mkey in mkeys:
                        metric_to_paths.setdefault(mkey, []).append(rel)
            except OSError:
                continue
    for mkey in list(metric_to_paths.keys()):
        metric_to_paths[mkey] = list(dict.fromkeys(metric_to_paths[mkey]))
    return metric_to_paths


# ══════════════════════════════════════════════════════════════
# Result JSON Generation
# ══════════════════════════════════════════════════════════════


def _coalesce_inheritance_from_diseases(var: Dict[str, Any]) -> str:
    """Top-level inheritance from annotator; else first non-empty diseases[].inheritance."""
    inh = (var.get("inheritance") or "").strip()
    if inh:
        return inh
    for d in var.get("diseases") or []:
        if isinstance(d, dict):
            x = (d.get("inheritance") or "").strip()
            if x:
                return x
    return ""


def _build_disease_groups_for_portal(
    variants: List[Dict[str, Any]],
    disease_variant_groups: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Portal Disease Summary: reviewData.disease_groups[disease] =
    { inheritance, carrier_frequency, severity, category, variants }.
    """
    vid_to_v = {v.get("variant_id"): v for v in variants if v.get("variant_id")}
    out: Dict[str, Any] = {}
    for g in disease_variant_groups:
        disease = g.get("disease") or ""
        if not disease:
            continue
        vobjs = [
            vid_to_v[vid]
            for vid in g.get("variant_ids") or []
            if vid in vid_to_v
        ]
        out[disease] = {
            "inheritance": g.get("inheritance") or "",
            "carrier_frequency": g.get("carrier_frequency") or "",
            "severity": g.get("severity") or "",
            "category": g.get("category") or "",
            "variants": vobjs,
        }
    return out


def generate_result_json(
    annotated_variants: List[Dict[str, Any]],
    acmg_results: List[Dict[str, Any]],
    qc_summary: Dict[str, Any],
    sample_name: str,
    order_id: str,
    disease_bed_info: Optional[Dict[str, Any]] = None,
    filter_summary: Optional[Dict[str, Any]] = None,
    parse_stats: Optional[Dict[str, Any]] = None,
    output_dir: str = "",
    metric_image_search_roots: Optional[List[str]] = None,
    analysis_dir: Optional[str] = None,
    dark_genes_extra_roots: Optional[List[str]] = None,
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
        filter_summary: 적용된 필터 요약 (선택)
        parse_stats: VCF 파싱 통계 (선택)
        output_dir: 출력 디렉토리
        metric_image_search_roots: QC 차트 이미지 탐색 루트 (artifact·layout output 등)
        analysis_dir: Pipeline analysis root (dark-gene Nextflow outdir); defaults to output_dir
        dark_genes_extra_roots: Extra dirs to scan for summary/*.txt when artifact analysis_dir
            has no pipeline publishDir output (e.g. Nextflow wrote under carrier_screening_layout_base).

    Returns:
        생성된 result.json 파일 경로
    """
    os.makedirs(output_dir, exist_ok=True)

    qc_summary["metric_images"] = discover_qc_metric_images(
        metric_image_search_roots or []
    )

    ad = (analysis_dir or output_dir or "").strip()
    od = (output_dir or "").strip()

    def _collect_dark(root: str) -> Dict[str, Any]:
        from .dark_genes import collect_dark_genes_from_analysis_dir

        return collect_dark_genes_from_analysis_dir(root, sample_name)

    search_roots: List[str] = []
    seen_real: set = set()

    def _add_dark_root(path: Optional[str]) -> None:
        if not path or not str(path).strip():
            return
        p = str(path).strip()
        if not os.path.isdir(p):
            return
        try:
            r = os.path.realpath(p)
        except OSError:
            return
        if r in seen_real:
            return
        seen_real.add(r)
        search_roots.append(os.path.abspath(p))

    try:
        _add_dark_root(ad)
        _add_dark_root(od)
        for extra in dark_genes_extra_roots or []:
            _add_dark_root(extra)
    except OSError:
        pass

    dark_genes_block: Dict[str, Any] = {
        "status": "not_found",
        "message": "No readable analysis/output directory for dark-gene summary",
        "searched_roots": search_roots,
    }
    if ad and not os.path.isdir(ad) and not search_roots:
        dark_genes_block["message"] = f"analysis_dir is not a directory: {ad}"

    try:
        if not search_roots:
            pass  # keep not_found above
        else:
            for root in search_roots:
                dark_genes_block = _collect_dark(root)
                if dark_genes_block.get("status") != "not_found":
                    break
            if dark_genes_block.get("status") == "not_found":
                dark_genes_block["searched_roots"] = search_roots
                dark_genes_block["message"] = (
                    "No dark_genes summary/detailed report under Nextflow outdir "
                    "(expected summary/*_summary_report.txt or *_summary_report.txt in outdir root). "
                    f"Tried: {search_roots}"
                )
    except Exception as e:
        logger.warning("[generate_result_json] dark_genes collection failed: %s", e)
        dark_genes_block = {"status": "error", "message": str(e)}

    try:
        prev_path = os.path.join(output_dir, "result.json")
        if os.path.isfile(prev_path):
            from .dark_genes import merge_dark_genes_reviews

            with open(prev_path, "r", encoding="utf-8") as rf:
                prev = json.load(rf)
            prev_dg = prev.get("dark_genes") if isinstance(prev, dict) else None
            if isinstance(prev_dg, dict):
                dark_genes_block = merge_dark_genes_reviews(dark_genes_block, prev_dg)
    except Exception as e:
        logger.warning("[generate_result_json] dark_genes section_reviews merge skipped: %s", e)

    # 변이 + ACMG 결합
    variants_for_review = []
    for i, (var, acmg) in enumerate(zip(annotated_variants, acmg_results)):
        entry = {
            "variant_id": f"VAR_{i+1:04d}",

            # 위치
            "chrom": var.get("chrom", ""),
            "pos": var.get("pos", ""),
            "ref": var.get("ref", ""),
            "alt": var.get("alt", ""),

            # 유전자/전사체
            "gene": var.get("gene", ""),
            "transcript": var.get("transcript", ""),
            "canonical_enst": var.get("canonical_enst", ""),
            "clinical_nm": var.get("clinical_nm", ""),
            "hgvsc": var.get("hgvsc", ""),
            "hgvsp": var.get("hgvsp", ""),
            "effect": var.get("effect", ""),

            # 샘플 메트릭스
            "dp": var.get("dp"),
            "ref_depth": var.get("ref_depth"),
            "alt_depth": var.get("alt_depth"),
            "vaf": var.get("vaf"),
            "gt": var.get("gt", ""),
            "zygosity": var.get("zygosity", ""),

            # gnomAD
            "gnomad_af": var.get("gnomad_af"),
            "gnomad_exomes_af": var.get("gnomad_exomes_af"),
            "gnomad_genomes_af": var.get("gnomad_genomes_af"),
            "gnomad_source": var.get("gnomad_source", ""),

            # ClinVar
            "clinvar_sig": var.get("clinvar_sig", ""),
            "clinvar_sig_primary": var.get("clinvar_sig_primary", ""),
            "clinvar_conflicting": var.get("clinvar_conflicting", False),
            "clinvar_conflict_detail": var.get("clinvar_conflict_detail", ""),
            "clinvar_revstat": var.get("clinvar_revstat", ""),
            "clinvar_stars": var.get("clinvar_stars", 0),
            "clinvar_dn": var.get("clinvar_dn", ""),
            "clinvar_variation_id": var.get("clinvar_variation_id", ""),

            # dbSNP
            "dbsnp_rsid": var.get("dbsnp_rsid", ""),
            "dbsnp_url": var.get("dbsnp_url", ""),

            # ClinGen
            "clingen_hi_score": var.get("clingen_hi_score"),
            "clingen_ts_score": var.get("clingen_ts_score"),
            "clingen_hi_desc": var.get("clingen_hi_desc", ""),
            "clingen_ts_desc": var.get("clingen_ts_desc", ""),
            "clingen_url": var.get("clingen_url", ""),

            # HGMD (라이선스 필요)
            "hgmd_class": var.get("hgmd_class", ""),
            "hgmd_disease": var.get("hgmd_disease", ""),
            "hgmd_pmid": var.get("hgmd_pmid", ""),

            # Curated Variant DB
            "curated_classification": var.get("curated_classification", ""),
            "curated_source": var.get("curated_source", ""),
            "curated_notes": var.get("curated_notes", ""),

            # Disease-Gene Mapping (inheritance duplicated at top level for reports / ClinVar-only rows)
            "diseases": var.get("diseases", []),
            "inheritance": _coalesce_inheritance_from_diseases(var),

            # HPO Phenotypes
            "hpo_phenotypes": var.get("hpo_phenotypes", []),

            # ACMG 분류
            "acmg_classification": acmg.get("final_classification", "VUS"),
            "acmg_criteria": acmg.get("final_criteria", []),
            "acmg_reasoning": acmg.get("final_reasoning", ""),
            "acmg_rule_based": acmg.get("rule_based", {}),
            "acmg_ai": acmg.get("ai"),
            "acmg_literature_enhanced": acmg.get("literature_enhanced", False),

            # 문헌 (PubMed) — pipeline이 literature_enabled 상태에서 AI 실행 시 채워짐
            # {total_found, from_cache, pmids, top_titles}
            "literature": acmg.get("literature"),

            # 리뷰어 입력 필드 (Portal에서 채움)
            "reviewer_classification": None,
            "reviewer_comment": None,
            "reviewer_confirmed": False,
            "include_in_report": False,
        }
        variants_for_review.append(entry)

    # 변이 통계
    variant_stats = _compute_variant_stats(variants_for_review)

    # IGV 스냅샷 경로 매핑
    igv_snapshots = _find_igv_snapshots(output_dir, variants_for_review)

    # 질환별 변이 그룹핑
    disease_variant_groups = _group_variants_by_disease(variants_for_review)
    disease_groups = _build_disease_groups_for_portal(variants_for_review, disease_variant_groups)

    # 유전자별 변이 그룹핑
    gene_variant_groups = _group_variants_by_gene(variants_for_review)

    result = {
        "version": "2.0",
        "type": "carrier_screening_result",
        "generated_at": now_kst_iso(),
        "order_id": order_id,
        "sample_name": sample_name,
        "status": "pending_review",

        # QC 요약
        "qc_summary": qc_summary,

        # 변이 통계
        "variant_stats": variant_stats,

        # 변이 목록 (리뷰 대상)
        "variants": variants_for_review,

        # 질환별 변이 그룹
        "disease_variant_groups": disease_variant_groups,
        # Portal Disease Summary (inheritance + panel metadata + variant objects)
        "disease_groups": disease_groups,

        # 유전자별 변이 그룹
        "gene_variant_groups": gene_variant_groups,

        # 질환 BED 정보
        "disease_panel": disease_bed_info or {},

        # 필터 요약
        "filter_summary": filter_summary or {},

        # 파싱 통계
        "parse_stats": parse_stats or {},

        # IGV 스냅샷 매핑
        "igv_snapshots": igv_snapshots,

        # Dark-gene / hard-to-sequence supplementary pipeline (Nextflow summary/)
        "dark_genes": dark_genes_block,

        # 메타데이터
        "metadata": {
            "pipeline_version": "carrier-screening-v2.0",
            "annotation_databases": {
                "clinvar": True,
                "gnomad": True,
                "snpeff": True,
                "clingen": True,
                "mane": True,
                "hpo": True,
                "hgmd": False,  # 라이선스 필요
                "curated_db": True,
                "disease_gene": True,
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
    genes = {}

    for v in variants:
        cls = v.get("acmg_classification", "VUS")
        classifications[cls] = classifications.get(cls, 0) + 1

        zyg = v.get("zygosity", "Unknown")
        zygosities[zyg] = zygosities.get(zyg, 0) + 1

        eff = v.get("effect", "Unknown")
        effects[eff] = effects.get(eff, 0) + 1

        gene = v.get("gene", "Unknown")
        genes[gene] = genes.get(gene, 0) + 1

    pathogenic_count = sum(
        classifications.get(c, 0)
        for c in ("Pathogenic", "Likely Pathogenic")
    )

    return {
        "total": total,
        "pathogenic_or_likely": pathogenic_count,
        "vus": classifications.get("VUS", 0),
        "benign_or_likely": sum(
            classifications.get(c, 0)
            for c in ("Benign", "Likely Benign")
        ),
        "by_classification": classifications,
        "by_zygosity": zygosities,
        "by_effect": effects,
        "by_gene": genes,
        "unique_genes": len(genes),
    }


def _group_variants_by_disease(variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """변이를 질환별로 그룹핑합니다 (질환별 inheritance·carrier_frequency 등은 첫 매칭 disease dict에서)."""
    disease_map: Dict[str, List[str]] = {}
    disease_meta: Dict[str, Dict[str, str]] = {}

    def _meta_from_disease_dict(d: Dict[str, Any]) -> Dict[str, str]:
        return {
            "inheritance": (d.get("inheritance") or "").strip(),
            "carrier_frequency": (d.get("carrier_frequency") or "").strip(),
            "severity": (d.get("severity") or "").strip(),
            "category": (d.get("category") or "").strip(),
        }

    for v in variants:
        diseases = v.get("diseases", [])
        clinvar_dn = v.get("clinvar_dn", "")
        variant_id = v.get("variant_id", "")

        if diseases:
            for d in diseases:
                name = (d.get("name") or d.get("disease_name") or "Unknown") if isinstance(d, dict) else str(d)
                disease_map.setdefault(name, []).append(variant_id)
                if name not in disease_meta and isinstance(d, dict):
                    disease_meta[name] = _meta_from_disease_dict(d)
        elif clinvar_dn and clinvar_dn not in (".", "not_provided", "not_specified"):
            disease_map.setdefault(clinvar_dn, []).append(variant_id)
            if clinvar_dn not in disease_meta:
                disease_meta[clinvar_dn] = {
                    "inheritance": (v.get("inheritance") or "").strip(),
                    "carrier_frequency": "",
                    "severity": "",
                    "category": "",
                }

    groups = []
    for disease, variant_ids in disease_map.items():
        meta = disease_meta.get(disease, {})
        groups.append({
            "disease": disease,
            "variant_ids": variant_ids,
            "variant_count": len(variant_ids),
            "inheritance": meta.get("inheritance", ""),
            "carrier_frequency": meta.get("carrier_frequency", ""),
            "severity": meta.get("severity", ""),
            "category": meta.get("category", ""),
        })

    groups.sort(key=lambda g: g["variant_count"], reverse=True)
    return groups


def _group_variants_by_gene(variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """변이를 유전자별로 그룹핑합니다."""
    gene_map: Dict[str, List[str]] = {}

    for v in variants:
        gene = v.get("gene", "Unknown")
        variant_id = v.get("variant_id", "")
        gene_map.setdefault(gene, []).append(variant_id)

    groups = []
    for gene, variant_ids in gene_map.items():
        groups.append({
            "gene": gene,
            "variant_ids": variant_ids,
            "variant_count": len(variant_ids),
        })

    groups.sort(key=lambda g: g["variant_count"], reverse=True)
    return groups


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
    "gene", "transcript", "clinical_nm", "hgvsc", "hgvsp", "effect",
    "zygosity", "gt", "dp", "vaf",
    "gnomad_af", "clinvar_sig_primary", "clinvar_stars",
    "clinvar_dn", "clinvar_conflicting", "clinvar_conflict_detail",
    "acmg_classification", "acmg_criteria",
    "dbsnp_rsid",
    "hgmd_class", "hgmd_disease",
    "curated_classification",
    "inheritance",
    "clingen_hi_score", "clingen_ts_score",
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
                row["acmg_criteria"] = ",".join(str(c) for c in criteria)
            # diseases를 문자열로 변환
            diseases = row.get("diseases", [])
            if isinstance(diseases, list):
                row["diseases"] = "; ".join(
                    (d.get("name") or d.get("disease_name") or "") for d in diseases if isinstance(d, dict)
                )
            writer.writerow(row)

    logger.info(f"Generated variants.tsv: {output_path}")
    return output_path
