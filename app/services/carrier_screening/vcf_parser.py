"""
VCF Parser Module

VCF 파일에서 변이를 파싱하고 기본 정보를 추출합니다.
phenotype_portal/main.py의 VCF 파싱 및 필터링 로직을 모듈화한 것입니다.

주요 기능:
    - VCF FORMAT 문제 필드 정리 (GQ, SB 등)
    - snpEff ANN/CSQ annotation 레이아웃 감지 및 파싱
    - 유전자, 전사체, HGVS, 기능적 영향 추출
    - 샘플 메트릭스 추출 (DP, AD, VAF, GT, Zygosity)
    - BED 기반 변이 필터링
    - HPO 기반 유전자 필터링
    - gnomAD max AF 필터링
    - ClinVar significance 필터링
    - ACMG classification 필터링
    - 통합 VCF 파싱 파이프라인 (parse_vcf_variants)
"""

import os
import re
import gzip
import logging
import asyncio
import subprocess
from typing import Dict, Any, List, Optional, Set, Tuple
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

# 표준 VEP CSQ 열 순서 (헤더 Format 파싱 실패 시 fall back)
_VEP_CSQ_DEFAULT_INDICES = (3, 6, 10, 11, 1)  # SYMBOL, Feature, HGVSc, HGVSp, Consequence

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

async def run_snpeff(
    input_vcf: str,
    output_vcf: str,
    snpeff_jar: str,
    snpeff_db: str = "GRCh38.86",
    data_dir: Optional[str] = None,
):
    """snpEff를 실행하여 VCF에 기능적 annotation을 추가합니다.

    async subprocess를 사용해 event loop를 블로킹하지 않습니다.
    동시에 여러 샘플의 annotation이 가능합니다.
    """
    if not os.path.exists(snpeff_jar):
        raise RuntimeError(f"snpEff jar not found at '{snpeff_jar}'")

    # snpEff는 snpEff.config 를 JAR와 같은 디렉터리(또는 cwd)에서 찾음
    snpeff_home = os.path.abspath(os.path.dirname(snpeff_jar))
    cfg = os.path.join(snpeff_home, "snpEff.config")
    if not os.path.isfile(cfg):
        raise RuntimeError(
            f"snpEff config missing: expected '{cfg}' next to the JAR. "
            "Copy snpEff.config from the snpEff core distribution into tools/snpEff/."
        )

    cmd = ["java", "-Xmx8g", "-jar", os.path.abspath(snpeff_jar)]
    if data_dir and str(data_dir).strip():
        d = os.path.abspath(str(data_dir).strip())
        os.makedirs(d, exist_ok=True)
        cmd.extend(["-dataDir", d])
    # 요약 HTML/텍스트는 cwd에 생성됨. Docker에서 JAR 디렉터리가 :ro 마운트면 실패하므로 끔.
    cmd.extend(["-noStats", snpeff_db, "-hgvs", input_vcf])
    logger.info(f"Running snpEff (cwd={snpeff_home}): {' '.join(cmd)}")

    with open(output_vcf, "w", encoding="utf-8") as out:
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            cwd=snpeff_home,
            stdout=out,
            stderr=asyncio.subprocess.PIPE,
        )
        _, stderr_bytes = await proc.communicate()

    if proc.returncode != 0:
        stderr_text = (stderr_bytes or b"").decode("utf-8", errors="replace")
        raise RuntimeError(f"snpEff error (exit {proc.returncode}):\n{stderr_text}")

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
        fields_part = ""
        idxfmt = desc.lower().find("format:")
        if idxfmt >= 0:
            tail = desc[idxfmt + len("format:") :].strip()
            if tail.startswith('"'):
                tail = tail[1:]
            fields_part = tail.split('"')[0].strip()
        fields = [f.strip() for f in fields_part.split("|")] if fields_part else []
        uf = [f.upper() for f in fields]

        def idx(names: Set[str]) -> Optional[int]:
            for n in names:
                if n in uf:
                    return uf.index(n)
            return None

        # SYMBOL must win over Gene (ENSG…): set iteration order is not stable across runtimes.
        gene_idx = idx({"SYMBOL"})
        if gene_idx is None:
            gene_idx = idx({"GENE"})
        tx_idx = idx({"FEATURE", "TRANSCRIPT"})
        c_idx = idx({"HGVSC", "HGVS_C"})
        p_idx = idx({"HGVSP", "HGVS_P"})
        eff_idx = idx({"CONSEQUENCE"})
        dg, dt, dc, dp, de = _VEP_CSQ_DEFAULT_INDICES
        if gene_idx is None:
            gene_idx = dg
        if tx_idx is None:
            tx_idx = dt
        if c_idx is None:
            c_idx = dc
        if p_idx is None:
            p_idx = dp
        if eff_idx is None:
            eff_idx = de
        return "CSQ", gene_idx, tx_idx, c_idx, p_idx, eff_idx

    return None, None, None, None, None, None


def _annotation_from_pipe_parts(
    parts: List[str],
    gi: Optional[int],
    ti: Optional[int],
    ci: Optional[int],
    pi: Optional[int],
    ei: Optional[int],
) -> Tuple[str, str, str, str, str]:
    """One CSQ/ANN pipe field list → gene, transcript, hgvsc, hgvsp, effect."""
    gene = transcript = hgvsc = hgvsp = effect = ""
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
    return gene, transcript, hgvsc, hgvsp, effect


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
        if tag == "CSQ":
            if isinstance(raw, (list, tuple)):
                csq_joined = ",".join(str(x) for x in raw)
            else:
                csq_joined = str(raw)
            picked = False
            for entry in csq_joined.split(","):
                entry = entry.strip()
                if not entry:
                    continue
                parts = entry.split("|")
                g, t, hc, hp, ef = _annotation_from_pipe_parts(parts, gi, ti, ci, pi, ei)
                if g:
                    gene, transcript, hgvsc, hgvsp, effect = g, t, hc, hp, ef
                    picked = True
                    break
            if not picked and csq_joined:
                first = csq_joined.split(",")[0].strip()
                if first:
                    gene, transcript, hgvsc, hgvsp, effect = _annotation_from_pipe_parts(
                        first.split("|"), gi, ti, ci, pi, ei
                    )
        else:
            if isinstance(raw, (list, tuple)):
                raw = raw[0]
            parts = str(raw).split("|")
            gene, transcript, hgvsc, hgvsp, effect = _annotation_from_pipe_parts(
                parts, gi, ti, ci, pi, ei
            )

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


def get_bed_gene_at_position(chrom: str, pos: int, bed_regions: Dict) -> str:
    """BED 영역에서 변이 위치의 유전자명을 반환합니다."""
    if not bed_regions:
        return ""
    regions = bed_regions.get(chrom, [])
    for start, end, name in regions:
        if start <= pos <= end:
            return name
    return ""


# ══════════════════════════════════════════════════════════════
# Filter Configuration (phenotype_portal 필터링 파라미터 통합)
# ══════════════════════════════════════════════════════════════

@dataclass
class VariantFilterConfig:
    """
    VCF 변이 필터링 설정.
    phenotype_portal/main.py의 /analyze 라우트 파라미터를 구조화한 것입니다.

    Attributes:
        hpo_genes: HPO 기반 유전자 집합 (비어있으면 필터 비활성)
        gene_filter_set: 수동 유전자 필터 집합 (비어있으면 필터 비활성)
        max_af: gnomAD AF 상한 (None이면 필터 비활성)
        clinvar_filter: ClinVar significance 필터 집합
            - "any": 모든 ClinVar 변이
            - "has": ClinVar에 등록된 변이
            - "plp": Pathogenic/Likely pathogenic
            - "vus": VUS
            - "blb": Benign/Likely benign
            - "conflicting": Conflicting
            - "not_conflicting": ClinVar에 있지만 conflict 아닌 것
        exclude_clinvar_conflicts: ClinVar conflict 변이 제외 여부
        acmg_filter: ACMG 분류 필터 집합 (비어있으면 필터 비활성)
        require_protein_altering: 단백질 변경 변이만 포함 여부
        backbone_bed_regions: backbone BED 영역 (비어있으면 필터 비활성)
        disease_bed_regions: disease BED 영역 (정보 추가용)
    """
    hpo_genes: Set[str] = field(default_factory=set)
    gene_filter_set: Set[str] = field(default_factory=set)
    max_af: Optional[float] = None
    clinvar_filter: Set[str] = field(default_factory=set)
    exclude_clinvar_conflicts: bool = False
    acmg_filter: Set[str] = field(default_factory=set)
    require_protein_altering: bool = True
    backbone_bed_regions: Dict = field(default_factory=dict)
    disease_bed_regions: Dict = field(default_factory=dict)


def apply_clinvar_filter(
    cv_primary: str,
    has_clinvar: bool,
    cv_conflict: bool,
    clinvar_filter: Set[str],
    exclude_conflicts: bool,
) -> bool:
    """
    ClinVar 필터를 적용합니다.
    phenotype_portal/main.py의 ClinVar 필터링 로직을 그대로 포팅한 것입니다.

    Returns:
        True: 변이를 포함, False: 변이를 제외
    """
    # conflict 제외 옵션
    if exclude_conflicts and cv_conflict:
        return False

    # ClinVar 필터가 없거나 "any"이면 모든 변이 통과
    if not clinvar_filter or "any" in clinvar_filter:
        return True

    matches = False
    if "has" in clinvar_filter and has_clinvar:
        matches = True
    if "plp" in clinvar_filter and cv_primary in ("Pathogenic", "Likely pathogenic"):
        matches = True
    if "vus" in clinvar_filter and cv_primary == "VUS":
        matches = True
    if "blb" in clinvar_filter and cv_primary in ("Benign", "Likely benign"):
        matches = True
    if "conflicting" in clinvar_filter and cv_conflict:
        matches = True
    if "not_conflicting" in clinvar_filter and has_clinvar and not cv_conflict:
        matches = True

    return matches


def _is_ensembl_gene_id(s: str) -> bool:
    t = (s or "").strip()
    return len(t) >= 12 and t.startswith("ENSG")


def _merge_ann_gene_from_csq_row(ann: Dict[str, Any], g2: str) -> None:
    """Fill or upgrade ann['gene'] using per-record CSQ (prefer HGNC symbol over ENSG)."""
    if not g2 or not str(g2).strip():
        return
    sym = str(g2).strip()
    cur = (ann.get("gene") or "").strip()
    if not cur:
        ann["gene"] = sym
        return
    if _is_ensembl_gene_id(cur) and not _is_ensembl_gene_id(sym):
        ann["gene"] = sym


def _lookup_vep_annotations(
    vep_annotations: Dict[str, Any],
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[Dict[str, Any]]:
    """
    Match pre-parsed VEP rows to pysam records.

    Keys from extract_vep_annotations_from_vcf use rec.chrom as in file (chr1 vs 1).
    """
    chroms: List[str] = [chrom]
    if str(chrom).startswith("chr") and len(str(chrom)) > 3:
        chroms.append(str(chrom)[3:])
    elif chrom and not str(chrom).startswith("chr"):
        chroms.append(f"chr{chrom}")

    def _ralt(x: str) -> List[str]:
        if not x:
            return [x]
        if x.startswith("<"):
            return [x]
        out = [x]
        u = x.upper()
        if u != x:
            out.append(u)
        return out

    keys: List[str] = []
    for c in chroms:
        for r in _ralt(str(ref)):
            for a in _ralt(str(alt)):
                keys.append(f"{c}:{pos}:{r}:{a}")
    for k in dict.fromkeys(keys):
        hit = vep_annotations.get(k)
        if hit:
            return hit
    return None


# ══════════════════════════════════════════════════════════════
# Integrated VCF Parsing Pipeline
# ══════════════════════════════════════════════════════════════

def parse_vcf_variants(
    vcf_path: str,
    annotator,
    filter_config: Optional[VariantFilterConfig] = None,
    acmg_classifier=None,
    vep_annotations: Optional[Dict[str, Any]] = None,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    """
    VCF 파일을 파싱하고 필터링/annotation을 수행하는 통합 파이프라인.
    phenotype_portal/main.py의 /analyze 라우트 로직을 모듈화한 것입니다.

    Args:
        vcf_path: 파싱할 VCF 파일 경로
        annotator: VariantAnnotator 인스턴스
        filter_config: 필터링 설정 (None이면 기본 설정 사용)
        acmg_classifier: ACMG 분류기 (None이면 rule-based만 사용)

    Returns:
        (annotated_variants, acmg_results, parse_stats)
        - annotated_variants: annotation된 변이 딕셔너리 리스트
        - acmg_results: ACMG 분류 결과 리스트
        - parse_stats: 파싱 통계 (total_records, parsed_genes, filtered_count, warnings)
    """
    import pysam

    if filter_config is None:
        filter_config = VariantFilterConfig()

    stats = {
        "total_records": 0,
        "parsed_gene_count": 0,
        "nonref_count": 0,
        "protein_altering_count": 0,
        "bed_filtered_count": 0,
        "gene_filtered_count": 0,
        "af_filtered_count": 0,
        "clinvar_filtered_count": 0,
        "acmg_filtered_count": 0,
        "final_count": 0,
        "warnings": [],
    }

    annotated_variants: List[Dict[str, Any]] = []
    acmg_results: List[Dict[str, Any]] = []

    try:
        vcf = pysam.VariantFile(vcf_path)
    except Exception as e:
        stats["warnings"].append(f"Failed to open VCF: {e}")
        return annotated_variants, acmg_results, stats

    # VEP annotation 모드 감지
    # vep_annotations가 제공된 경우 VEP CSQ 기반 annotation 경로를 사용합니다.
    # 제공되지 않은 경우 기존 snpEff ANN/CSQ 파싱 경로를 사용합니다.
    use_vep_path = bool(vep_annotations)
    layout_tag, layout_gi, layout_ti, layout_ci, layout_pi, layout_ei = get_annotation_layout(vcf)
    if use_vep_path:
        logger.info("VEP annotation mode: using pre-parsed CSQ data (%d variants)", len(vep_annotations))
        # Still use header CSQ/ANN layout for the *per-record* pass so gene / effect are populated
        # before HPO and gene_filter_set. Previously tag was cleared here, leaving gene="" and dropping
        # every record when those filters were configured (dark-gene / phenotype orders).
        if layout_tag:
            tag, gi, ti, ci, pi, ei = (
                layout_tag,
                layout_gi,
                layout_ti,
                layout_ci,
                layout_pi,
                layout_ei,
            )
        else:
            tag, gi, ti, ci, pi, ei = None, None, None, None, None, None
    else:
        tag, gi, ti, ci, pi, ei = (
            layout_tag,
            layout_gi,
            layout_ti,
            layout_ci,
            layout_pi,
            layout_ei,
        )
        if tag is None:
            stats["warnings"].append("No ANN/CSQ detected; c./p. may be missing.")

    gene_interval_lookup = None
    if annotator and hasattr(annotator, 'gene_intervals') and annotator.gene_intervals:
        gene_interval_lookup = annotator.gene_intervals.lookup

    try:
        for rec in vcf:
            stats["total_records"] += 1

            # 1. non-reference 필터
            if not is_nonref(rec):
                continue
            stats["nonref_count"] += 1

            # 2. 유전자/전사체/HGVS 추출
            gene, transcript, hgvsc, hgvsp, effect = extract_variant_info(
                rec, tag, gi, ti, ci, pi, ei,
                gene_interval_lookup=gene_interval_lookup,
            )
            if gene:
                stats["parsed_gene_count"] += 1

            # 3. 단백질 변경 필터
            if filter_config.require_protein_altering and not is_protein_altering(effect):
                continue
            stats["protein_altering_count"] += 1

            # 4. BED 영역 필터 (backbone)
            if filter_config.backbone_bed_regions:
                if not variant_in_bed(rec.chrom, rec.pos, filter_config.backbone_bed_regions):
                    stats["bed_filtered_count"] += 1
                    continue

            # 5. HPO 유전자 필터
            if filter_config.hpo_genes and gene not in filter_config.hpo_genes:
                stats["gene_filtered_count"] += 1
                continue

            # 6. 수동 유전자 필터
            if filter_config.gene_filter_set and gene.upper() not in filter_config.gene_filter_set:
                stats["gene_filtered_count"] += 1
                continue

            # 각 ALT allele에 대해 처리
            alts = list(rec.alts or [])
            if not alts:
                continue

            for alt in alts:
                alt_str = str(alt)

                # 샘플 메트릭스 (rec.id 를 함께 전달해 rsID 폴백에 활용)
                sample_metrics = get_sample_metrics(rec, alt_str)
                rec_id = rec.id if rec.id and rec.id != "." else ""
                if rec_id:
                    sample_metrics["rec_id"] = rec_id

                # 통합 annotation
                # VEP annotation이 제공된 경우 VEP 기반 경로 사용
                vep_data = None
                if use_vep_path and vep_annotations:
                    vep_data = _lookup_vep_annotations(
                        vep_annotations, rec.chrom, rec.pos, rec.ref, alt_str
                    )
                if vep_data and hasattr(annotator, 'annotate_with_vep'):
                    ann = annotator.annotate_with_vep(
                        chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=alt_str,
                        vep_data=vep_data,
                        sample_metrics=sample_metrics,
                    )
                    # VEP에서 유전자/전사체 정보 보완 (BED 필터용)
                    if not gene and ann.get("gene"):
                        gene = ann["gene"]
                    # Pre-parsed CSQ row can lack SYMBOL/HGVSc (wrong transcript picked, header drift).
                    # Merge missing fields from per-record CSQ/ANN so review table is not blank.
                    if use_vep_path and layout_tag:
                        g2, t2, hc2, hp2, ef2 = extract_variant_info(
                            rec,
                            layout_tag,
                            layout_gi,
                            layout_ti,
                            layout_ci,
                            layout_pi,
                            layout_ei,
                            gene_interval_lookup=gene_interval_lookup,
                        )
                        _merge_ann_gene_from_csq_row(ann, g2)
                        if t2 and not (ann.get("transcript") or "").strip():
                            ann["transcript"] = t2
                        if hc2 and not (ann.get("hgvsc") or "").strip():
                            ann["hgvsc"] = hc2
                        if hp2 and not (ann.get("hgvsp") or "").strip():
                            ann["hgvsp"] = hp2
                        if ef2 and not (ann.get("effect") or "").strip():
                            ann["effect"] = ef2
                    if ann.get("gnomad_af") is None and hasattr(annotator, "gnomad") and annotator.gnomad:
                        try:
                            gm = annotator.gnomad.lookup(rec.chrom, rec.pos, rec.ref, alt_str)
                            if gm and gm.get("af") is not None:
                                ann["gnomad_af"] = gm.get("af")
                                if ann.get("gnomad_exomes_af") is None:
                                    ann["gnomad_exomes_af"] = gm.get("exomes_af")
                                if ann.get("gnomad_genomes_af") is None:
                                    ann["gnomad_genomes_af"] = gm.get("genomes_af")
                                if not (ann.get("gnomad_source") or "").strip():
                                    ann["gnomad_source"] = "local_vcf"
                        except Exception:
                            pass
                else:
                    g2, t2, hc2, hp2, ef2 = gene, transcript, hgvsc, hgvsp, effect
                    # VEP pre-parse missed this row (key mismatch) or non-VEP VCF: parse ANN/CSQ inline
                    if use_vep_path and layout_tag:
                        g2, t2, hc2, hp2, ef2 = extract_variant_info(
                            rec,
                            layout_tag,
                            layout_gi,
                            layout_ti,
                            layout_ci,
                            layout_pi,
                            layout_ei,
                            gene_interval_lookup=gene_interval_lookup,
                        )
                    ann = annotator.annotate(
                        chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=alt_str,
                        gene=g2, transcript=t2,
                        hgvsc=hc2, hgvsp=hp2, effect=ef2,
                        sample_metrics=sample_metrics,
                    )

                # 7. gnomAD max AF 필터
                gnomad_af = ann.get("gnomad_af")
                if filter_config.max_af is not None and gnomad_af is not None:
                    if gnomad_af > filter_config.max_af:
                        stats["af_filtered_count"] += 1
                        continue

                # 8. ClinVar 필터
                cv_primary = ann.get("clinvar_sig_primary", "")
                has_clinvar = bool(
                    ann.get("clinvar_sig") or
                    ann.get("clinvar_variation_id") or
                    ann.get("dbsnp_rsid")
                )
                cv_conflict = bool(ann.get("clinvar_conflicting", False))

                if filter_config.clinvar_filter:
                    if not apply_clinvar_filter(
                        cv_primary, has_clinvar, cv_conflict,
                        filter_config.clinvar_filter,
                        filter_config.exclude_clinvar_conflicts,
                    ):
                        stats["clinvar_filtered_count"] += 1
                        continue
                elif filter_config.exclude_clinvar_conflicts and cv_conflict:
                    stats["clinvar_filtered_count"] += 1
                    continue

                # 9. ACMG 분류
                acmg = {"final_classification": "VUS", "final_criteria": [], "final_reasoning": ""}
                if acmg_classifier:
                    try:
                        ctx = {
                            "af": gnomad_af,
                            "clinvar": cv_primary or ann.get("clinvar_sig", ""),
                            "effect": ann.get("effect") or effect,
                            "chrom": rec.chrom,
                            "pos": rec.pos,
                            "ref": rec.ref,
                            "alt": alt_str,
                        }
                        result = acmg_classifier.classify(
                            ann.get("gene") or gene,
                            ann.get("hgvsc") or hgvsc,
                            context=ctx,
                            lite_mode=True, do_local_lookup=False,
                        )
                        acmg = {
                            "final_classification": getattr(result, "acmg_classification", "VUS"),
                            "final_criteria": getattr(result, "acmg_evidence", []),
                            "final_reasoning": "",
                        }
                    except Exception as e:
                        logger.warning(f"ACMG classifier error for {gene}: {e}")

                # 10. ACMG 필터
                if filter_config.acmg_filter:
                    acmg_class = acmg.get("final_classification", "VUS")
                    if acmg_class not in filter_config.acmg_filter:
                        stats["acmg_filtered_count"] += 1
                        continue

                # Disease BED 정보 추가
                if filter_config.disease_bed_regions:
                    bed_gene = get_bed_gene_at_position(
                        rec.chrom, rec.pos, filter_config.disease_bed_regions
                    )
                    if bed_gene:
                        ann["disease_bed_gene"] = bed_gene

                annotated_variants.append(ann)
                acmg_results.append(acmg)

    except OSError as e:
        stats["warnings"].append(f"Error parsing VCF: {e}")

    vcf.close()

    stats["final_count"] = len(annotated_variants)
    logger.info(
        f"VCF parsing complete: {stats['total_records']} records → "
        f"{stats['final_count']} variants after filtering"
    )

    return annotated_variants, acmg_results, stats
