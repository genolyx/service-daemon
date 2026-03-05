"""
VCF Parser Module

VCF 파일에서 변이를 파싱하고 기본 정보를 추출합니다.
phenotype_portal/main.py의 VCF 파싱 로직을 모듈화한 것입니다.

주요 기능:
    - VCF FORMAT 문제 필드 정리 (GQ, SB 등)
    - snpEff ANN/CSQ annotation 레이아웃 감지 및 파싱
    - 유전자, 전사체, HGVS, 기능적 영향 추출
    - 샘플 메트릭스 추출 (DP, AD, VAF, GT, Zygosity)
    - BED 기반 변이 필터링
"""

import os
import re
import gzip
import logging
import subprocess
from typing import Dict, Any, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# ─── VCF FORMAT 문제 필드 정리 ─────────────────────────────

TROUBLE_FORMATS = {"GQ", "SB"}


def clean_vcf_remove_formats(input_vcf: str, output_vcf: str):
    """
    pysam 파싱 에러를 유발하는 FORMAT 필드(GQ, SB 등)를 제거합니다.
    '.' 값이 Integer 타입 필드에 들어가는 경우 pysam이 크래시하는 문제를 해결합니다.
    """
    opener_in = gzip.open if input_vcf.endswith(".gz") else open
    with opener_in(input_vcf, "rt", encoding="utf-8", errors="replace") as fin, \
            open(output_vcf, "w", encoding="utf-8") as fout:
        for line in fin:
            line = line.rstrip("\n")

            if line.startswith("##FORMAT=<ID="):
                if any(f"ID={fid}" in line for fid in TROUBLE_FORMATS):
                    continue
                fout.write(line + "\n")
                continue

            if line.startswith("#"):
                fout.write(line + "\n")
                continue

            parts = line.split("\t")
            if len(parts) <= 8:
                fout.write(line + "\n")
                continue

            fmt_fields = parts[8].split(":")
            trouble_indices = [i for i, name in enumerate(fmt_fields) if name in TROUBLE_FORMATS]
            if not trouble_indices:
                fout.write(line + "\n")
                continue

            keep_indices = [i for i in range(len(fmt_fields)) if i not in trouble_indices]
            parts[8] = ":".join(fmt_fields[i] for i in keep_indices) if keep_indices else "."

            for i in range(9, len(parts)):
                sample = parts[i]
                if sample in (".", ""):
                    continue
                sf = sample.split(":")
                if len(sf) != len(fmt_fields):
                    continue
                parts[i] = ":".join(sf[j] for j in keep_indices) if keep_indices else "."

            fout.write("\t".join(parts) + "\n")


# ─── snpEff Annotation ────────────────────────────────────

def run_snpeff(input_vcf: str, output_vcf: str, snpeff_jar: str, snpeff_db: str = "GRCh38.86"):
    """snpEff를 실행하여 VCF에 기능적 annotation을 추가합니다."""
    if not os.path.exists(snpeff_jar):
        raise RuntimeError(f"snpEff jar not found at '{snpeff_jar}'")

    cmd = ["java", "-Xmx4g", "-jar", snpeff_jar, snpeff_db, "-hgvs", input_vcf]
    logger.info(f"Running snpEff: {' '.join(cmd)}")

    with open(output_vcf, "w", encoding="utf-8") as out:
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)

    if res.returncode != 0:
        raise RuntimeError(f"snpEff error (exit {res.returncode}):\n{res.stderr}")

    logger.info(f"snpEff annotation complete: {output_vcf}")


# ─── Annotation Layout Detection (ANN/CSQ) ────────────────

def get_annotation_layout(vcf) -> Tuple:
    """
    VCF 헤더에서 ANN 또는 CSQ 필드의 레이아웃을 감지합니다.

    Returns:
        (tag, gene_idx, transcript_idx, hgvsc_idx, hgvsp_idx, effect_idx)
        tag가 None이면 annotation이 없는 VCF입니다.
    """
    if hasattr(vcf, 'header'):
        header_info = vcf.header.info
    else:
        return None, None, None, None, None, None

    if "ANN" in header_info:
        desc = header_info["ANN"].description or ""
        fields_part = desc.split("Format:")[-1].strip().strip('"')
        fields = [f.strip() for f in fields_part.split("|")] if fields_part else []
        uf = [f.upper() for f in fields]

        def idx(names: Set[str]) -> Optional[int]:
            for n in names:
                if n in uf:
                    return uf.index(n)
            return None

        gene_idx = idx({"GENE_NAME", "GENE"})
        tx_idx = idx({"FEATURE_ID", "TRANSCRIPT"})
        c_idx = idx({"HGVSC", "HGVS.C"})
        p_idx = idx({"HGVSP", "HGVS.P"})
        eff_idx = idx({"ANNOTATION"})

        # 기본 snpEff 인덱스 폴백
        if not fields or gene_idx is None or tx_idx is None or c_idx is None or p_idx is None or eff_idx is None:
            gene_idx = 3
            eff_idx = 1
            tx_idx = 6
            c_idx = 9
            p_idx = 10

        return "ANN", gene_idx, tx_idx, c_idx, p_idx, eff_idx

    if "CSQ" in header_info:
        desc = header_info["CSQ"].description or ""
        fields_part = desc.split("Format:")[-1].strip().strip('"')
        fields = [f.strip() for f in fields_part.split("|")] if fields_part else []
        uf = [f.upper() for f in fields]

        def idx(names: Set[str]) -> Optional[int]:
            for n in names:
                if n in uf:
                    return uf.index(n)
            return None

        gene_idx = idx({"SYMBOL", "GENE"})
        tx_idx = idx({"FEATURE", "TRANSCRIPT"})
        c_idx = idx({"HGVSC", "HGVS_C"})
        p_idx = idx({"HGVSP", "HGVS_P"})
        eff_idx = idx({"CONSEQUENCE"})
        return "CSQ", gene_idx, tx_idx, c_idx, p_idx, eff_idx

    return None, None, None, None, None, None


def extract_variant_info(rec, tag, gi, ti, ci, pi, ei, gene_interval_lookup=None) -> Tuple[str, str, str, str, str]:
    """
    VCF 레코드에서 유전자, 전사체, HGVS, 기능적 영향을 추출합니다.

    Args:
        rec: pysam VariantRecord
        tag, gi, ti, ci, pi, ei: annotation 레이아웃 인덱스
        gene_interval_lookup: 위치 기반 유전자 조회 함수 (선택)

    Returns:
        (gene, transcript, hgvsc, hgvsp, effect)
    """
    gene = transcript = hgvsc = hgvsp = effect = ""

    if tag and tag in rec.info:
        raw = rec.info[tag]
        if isinstance(raw, (list, tuple)):
            raw = raw[0]
        parts = str(raw).split("|")

        if gi is not None and gi < len(parts):
            g = parts[gi].strip()
            if g and g != ".":
                gene = g

        if ti is not None and ti < len(parts):
            t = parts[ti].strip()
            if t and t != ".":
                transcript = t

        if ci is not None and ci < len(parts):
            cv = parts[ci].strip()
            if cv and cv != ".":
                if ":" in cv:
                    tx, cval = cv.split(":", 1)
                    if not transcript:
                        transcript = tx.strip()
                    hgvsc = cval.strip()
                else:
                    hgvsc = cv

        if pi is not None and pi < len(parts):
            pv = parts[pi].strip()
            if pv and pv != ".":
                if ":" in pv:
                    hgvsp = pv.split(":", 1)[1].strip()
                else:
                    hgvsp = pv

        if ei is not None and ei < len(parts):
            ev = parts[ei].strip()
            if ev and ev != ".":
                effect = ev

    # 유전자를 찾지 못한 경우 위치 기반 조회
    if not gene and gene_interval_lookup:
        gene = gene_interval_lookup(rec.chrom, rec.pos) or ""

    return gene, transcript, hgvsc, hgvsp, effect


# ─── Sample Metrics Extraction ─────────────────────────────

def _safe_int(x) -> Optional[int]:
    try:
        if x is None:
            return None
        if isinstance(x, bool):
            return None
        return int(x)
    except Exception:
        try:
            return int(float(str(x)))
        except Exception:
            return None


def _safe_float(x) -> Optional[float]:
    try:
        if x is None:
            return None
        if isinstance(x, bool):
            return None
        return float(x)
    except Exception:
        try:
            return float(str(x))
        except Exception:
            return None


def _alts_as_str_list(rec) -> List[str]:
    return [str(a) for a in (rec.alts or [])]


def _pick_first_sample(rec):
    if not rec.samples:
        return None, None
    name = next(iter(rec.samples.keys()))
    return name, rec.samples[name]


def _zygosity_from_gt(gt: Tuple[Optional[int], ...]) -> Tuple[str, str]:
    if not gt:
        return "", ""
    gt_str = "/".join("." if a is None else str(a) for a in gt)
    if len(gt) < 2:
        return "Other", gt_str
    a0, a1 = gt[0], gt[1]
    if a0 is None or a1 is None:
        return "Other", gt_str
    if (a0 == 0 and a1 == 1) or (a0 == 1 and a1 == 0):
        return "Het", gt_str
    if a0 == 1 and a1 == 1:
        return "Hom", gt_str
    if a0 == 0 and a1 == 0:
        return "Ref", gt_str
    return "Other", gt_str


def get_sample_metrics(rec, alt_allele: str) -> Dict[str, Any]:
    """
    VCF 레코드에서 샘플 수준의 메트릭스를 추출합니다.

    Returns:
        {sample_name, dp, ref_depth, alt_depth, vaf, gt, zygosity, metrics_source}
    """
    dp = None
    ref_depth = None
    alt_depth = None
    vaf = None
    gt_str = ""
    zygosity = ""
    source = ""

    alts = _alts_as_str_list(rec)
    alt_index = 0
    if alts:
        try:
            alt_index = alts.index(str(alt_allele))
        except ValueError:
            alt_index = 0

    sample_name, s = _pick_first_sample(rec)

    if s is not None:
        if "DP" in s and s["DP"] is not None:
            dp = _safe_int(s["DP"])
            if dp is not None:
                source = source or "FORMAT:DP"

        gt = s.get("GT")
        if gt:
            zygosity, gt_str = _zygosity_from_gt(gt)

        if "AD" in s and s["AD"] is not None:
            try:
                ad = list(s["AD"])
                if len(ad) >= 2:
                    rd = _safe_int(ad[0])
                    ai = 1 + alt_index
                    ad_alt = _safe_int(ad[ai]) if ai < len(ad) else None
                    if rd is not None:
                        ref_depth = rd
                    if ad_alt is not None:
                        alt_depth = ad_alt
                    if ref_depth is not None and alt_depth is not None:
                        denom = ref_depth + alt_depth
                        if denom > 0:
                            vaf = alt_depth / denom
                            source = source or "FORMAT:AD"
            except Exception:
                pass

    # INFO:DP 폴백
    if dp is None and "DP" in rec.info and rec.info["DP"] is not None:
        val = rec.info["DP"]
        if isinstance(val, (list, tuple)):
            val = val[0] if val else None
        dp = _safe_int(val)
        if dp is not None:
            source = source or "INFO:DP"

    if vaf is None and ref_depth is not None and alt_depth is not None:
        denom = ref_depth + alt_depth
        if denom > 0:
            vaf = alt_depth / denom

    if dp is None and ref_depth is not None and alt_depth is not None:
        dp = ref_depth + alt_depth

    if zygosity == "Ref":
        zygosity = "Other"

    return {
        "sample_name": sample_name or "",
        "dp": dp,
        "ref_depth": ref_depth,
        "alt_depth": alt_depth,
        "vaf": round(vaf, 4) if vaf is not None else None,
        "gt": gt_str,
        "zygosity": zygosity,
        "metrics_source": source,
    }


# ─── Variant Filters ──────────────────────────────────────

PROTEIN_ALTERING = {
    "missense_variant", "stop_gained", "stop_lost",
    "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "start_lost",
    "inframe_insertion", "inframe_deletion",
    "protein_altering_variant",
}


def is_protein_altering(effect: str) -> bool:
    """변이가 단백질에 영향을 미치는지 확인합니다."""
    if not effect:
        return True  # annotation이 없으면 포함
    parts = re.split(r"[&,]", effect)
    return any(p.strip() in PROTEIN_ALTERING for p in parts)


def is_nonref(rec) -> bool:
    """샘플이 non-reference genotype인지 확인합니다."""
    if not rec.samples:
        return True
    for s in rec.samples.values():
        gt = s.get("GT")
        if gt and any(a != 0 for a in gt if a is not None):
            return True
    return False


# ─── BED Region Filter ────────────────────────────────────

def load_bed_regions(bed_path: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    BED 파일을 로드하여 {chrom: [(start, end, name), ...]} 형태로 반환합니다.
    """
    regions = {}
    if not bed_path or not os.path.exists(bed_path):
        logger.warning(f"BED file not found: {bed_path}")
        return regions

    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3].strip() if len(parts) > 3 else ""
            regions.setdefault(chrom, []).append((start, end, name))

    total = sum(len(v) for v in regions.values())
    logger.info(f"Loaded {total} regions from BED: {bed_path}")
    return regions


def variant_in_bed(chrom: str, pos: int, bed_regions: Dict) -> bool:
    """변이가 BED 영역 내에 있는지 확인합니다."""
    if not bed_regions:
        return True  # BED가 없으면 모든 변이 통과

    regions = bed_regions.get(chrom, [])
    for start, end, _ in regions:
        if start <= pos <= end:
            return True
    return False
