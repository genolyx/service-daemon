"""
vep_parser.py — VEP (Variant Effect Predictor) CSQ INFO 필드 파서

Nextflow 파이프라인에서 VEP annotation이 완료된 VCF를 service-daemon이 수신할 때,
INFO 필드의 CSQ 태그를 파싱하여 구조화된 딕셔너리로 변환합니다.

VEP CSQ 필드 형식:
    ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from
    Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|
    BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
    Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|
    MANE_SELECT|MANE_PLUS_CLINICAL|CANONICAL|SIFT|PolyPhen|AF|gnomADe_AF|gnomADg_AF|
    CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|
    TRANSCRIPTION_FACTORS">

이 모듈은 VEP가 없는 환경(skip_vep=true, snpEff 사용)에서도 graceful fallback을
제공하므로 기존 snpEff 기반 파이프라인과 완전히 호환됩니다.
"""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# VEP CSQ 필드 순서 (VEP --pick 옵션 기준 표준 출력 필드)
# 실제 VCF 헤더의 CSQ Description에서 파싱하는 것이 가장 정확하지만,
# 표준 필드 목록을 fallback으로 제공합니다.
_DEFAULT_CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
    "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON",
    "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
    "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND",
    "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "MANE_SELECT", "MANE_PLUS_CLINICAL",
    "CANONICAL", "SIFT", "PolyPhen", "AF", "gnomADe_AF", "gnomADg_AF",
    "CLIN_SIG", "SOMATIC", "PHENO",
]


def parse_csq_header(vcf_header_str: str) -> Optional[List[str]]:
    """
    VCF 헤더 문자열에서 CSQ 필드 목록을 추출합니다.

    Args:
        vcf_header_str: VCF 헤더 전체 문자열 (##INFO=<ID=CSQ,...> 라인 포함)

    Returns:
        CSQ 필드 이름 목록 또는 None (CSQ 헤더가 없는 경우)
    """
    # ##INFO=<ID=CSQ,...,Description="... Format: field1|field2|...">
    match = re.search(
        r'##INFO=<ID=CSQ[^>]*Description="[^"]*Format:\s*([^"]+)"',
        vcf_header_str,
        re.IGNORECASE,
    )
    if not match:
        return None
    fields_str = match.group(1).strip()
    return [f.strip() for f in fields_str.split("|")]


def parse_csq_record(
    csq_value: str,
    csq_fields: List[str],
    alt_allele: Optional[str] = None,
    prefer_mane: bool = True,
) -> Optional[Dict[str, Any]]:
    """
    단일 변이의 CSQ INFO 값을 파싱하여 가장 적합한 transcript annotation을 반환합니다.

    VEP --pick 옵션을 사용하면 이미 하나의 transcript만 선택되지만,
    여러 transcript가 있을 경우 MANE_SELECT > CANONICAL 순으로 우선순위를 적용합니다.

    Args:
        csq_value:   VCF INFO 필드의 CSQ 값 (쉼표로 구분된 여러 transcript 포함 가능)
        csq_fields:  CSQ 필드 이름 목록 (parse_csq_header로 추출)
        alt_allele:  현재 처리 중인 ALT allele (여러 ALT가 있을 때 필터링용)
        prefer_mane: True면 MANE Select transcript를 우선 선택

    Returns:
        파싱된 annotation 딕셔너리 또는 None
    """
    if not csq_value or not csq_fields:
        return None

    transcripts = []
    for entry in csq_value.split(","):
        parts = entry.split("|")
        # 필드 수가 맞지 않으면 빈 문자열로 패딩
        if len(parts) < len(csq_fields):
            parts.extend([""] * (len(csq_fields) - len(parts)))
        record = dict(zip(csq_fields, parts))
        transcripts.append(record)

    if not transcripts:
        return None

    # ALT allele 필터링 (multi-allelic 변이 처리)
    if alt_allele:
        filtered = [t for t in transcripts if t.get("Allele", "") == alt_allele]
        if filtered:
            transcripts = filtered

    # MANE Select 우선 선택
    if prefer_mane:
        mane_transcripts = [
            t for t in transcripts
            if t.get("MANE_SELECT", "") and t["MANE_SELECT"] != ""
        ]
        if mane_transcripts:
            transcripts = mane_transcripts
        else:
            # CANONICAL transcript 차선 선택
            canonical = [t for t in transcripts if t.get("CANONICAL", "") == "YES"]
            if canonical:
                transcripts = canonical

    # 가장 심각한 consequence를 가진 transcript 선택
    impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    transcripts.sort(key=lambda t: impact_order.get(t.get("IMPACT", "MODIFIER"), 3))
    best = transcripts[0]

    def _sym_ok(t: Dict[str, str]) -> bool:
        s = (t.get("SYMBOL") or "").strip()
        return bool(s and s != ".")

    # Top hit can be intergenic / regulatory with empty SYMBOL; prefer a gene-bearing CSQ row.
    if transcripts and not _sym_ok(best):
        with_sym = [t for t in transcripts if _sym_ok(t)]
        if with_sym:
            with_sym.sort(key=lambda t: impact_order.get(t.get("IMPACT", "MODIFIER"), 3))
            best = with_sym[0]

    return _normalize_csq_record(best)


def _normalize_csq_record(record: Dict[str, str]) -> Dict[str, Any]:
    """
    CSQ 레코드를 service-daemon의 annotator.py가 기대하는 형식으로 정규화합니다.

    annotator.py의 annotate() 메서드 반환값과 동일한 키 이름을 사용하여
    기존 vcf_parser.py, review.py 등과의 호환성을 유지합니다.
    """

    def _float_or_none(val: str) -> Optional[float]:
        try:
            return float(val) if val and val not in (".", "") else None
        except (ValueError, TypeError):
            return None

    def _clean(val: str) -> str:
        return val if val and val != "." else ""

    # gnomAD AF: exomes와 genomes 중 더 큰 값(max)을 최종 AF로 사용
    # 기존 service-daemon 로컬 gnomAD 조회 방식(GnomADAnnotator.lookup)과 동일한 기준
    # → AF 필터(max_af) 적용 시 더 보수적으로 동작하여 false positive를 최소화
    gnomad_exomes_af = _float_or_none(record.get("gnomADe_AF", ""))
    gnomad_genomes_af = _float_or_none(record.get("gnomADg_AF", ""))
    _af_vals = [x for x in (gnomad_exomes_af, gnomad_genomes_af) if x is not None]
    gnomad_af = max(_af_vals) if _af_vals else None

    # HGVSc: "ENST00000123456.1:c.123A>G" → "c.123A>G"
    hgvsc_raw = _clean(record.get("HGVSc", ""))
    hgvsc = hgvsc_raw.split(":")[-1] if ":" in hgvsc_raw else hgvsc_raw

    # HGVSp: "ENSP00000123456.1:p.Arg41Gln" → "p.Arg41Gln"
    hgvsp_raw = _clean(record.get("HGVSp", ""))
    hgvsp = hgvsp_raw.split(":")[-1] if ":" in hgvsp_raw else hgvsp_raw
    # VEP은 아미노산을 3글자 코드로 표기 (p.Arg41Gln)
    # 필요 시 1글자 코드로 변환 가능하지만 임상 보고서에서는 3글자 권장

    # Transcript (Feature): ENST → NM 변환은 MANE 매핑에서 처리
    transcript = _clean(record.get("Feature", ""))

    # dbSNP rsID: Existing_variation 필드에서 추출 (rs로 시작하는 항목)
    existing = _clean(record.get("Existing_variation", ""))
    dbsnp_rsid = ""
    if existing:
        for item in existing.split("&"):
            if item.startswith("rs"):
                dbsnp_rsid = item
                break

    # SIFT 파싱: "tolerated(0.23)" → score=0.23, pred="tolerated"
    sift_raw = _clean(record.get("SIFT", ""))
    sift_score, sift_pred = _parse_score_pred(sift_raw)

    # PolyPhen 파싱: "probably_damaging(0.987)" → score=0.987, pred="probably_damaging"
    polyphen_raw = _clean(record.get("PolyPhen", ""))
    polyphen_score, polyphen_pred = _parse_score_pred(polyphen_raw)

    return {
        # 유전자/전사체 정보 (annotator.py 호환 키)
        # SYMBOL empty: fall back to VEP Gene (Ensembl gene id, e.g. ENSG…) for review visibility
        "gene": _clean(record.get("SYMBOL", "")) or _clean(record.get("Gene", "")),
        "transcript": transcript,
        "hgvsc": hgvsc,
        "hgvsp": hgvsp,
        "effect": _clean(record.get("Consequence", "")),
        "impact": _clean(record.get("IMPACT", "")),
        "biotype": _clean(record.get("BIOTYPE", "")),
        "exon": _clean(record.get("EXON", "")),
        "intron": _clean(record.get("INTRON", "")),
        # MANE / Canonical 정보
        "mane_select": _clean(record.get("MANE_SELECT", "")),
        "mane_plus_clinical": _clean(record.get("MANE_PLUS_CLINICAL", "")),
        "is_canonical": record.get("CANONICAL", "") == "YES",
        # 집단 빈도 (gnomAD)
        "gnomad_af": gnomad_af,
        "gnomad_exomes_af": gnomad_exomes_af,
        "gnomad_genomes_af": gnomad_genomes_af,
        "gnomad_source": "vep_csq",
        # dbSNP
        "dbsnp_rsid": dbsnp_rsid,
        # 기능 예측 점수
        "sift_score": sift_score,
        "sift_pred": sift_pred,
        "polyphen_score": polyphen_score,
        "polyphen_pred": polyphen_pred,
        # ClinVar (VEP에서 제공하는 CLIN_SIG; 정밀 매칭은 daemon annotator에서 수행)
        "vep_clin_sig": _clean(record.get("CLIN_SIG", "")),
        # 원본 CSQ 레코드 보존 (디버깅용)
        "_csq_raw": record,
    }


def _parse_score_pred(raw: str):
    """
    "probably_damaging(0.987)" 형식의 VEP 예측값을 (score, prediction) 튜플로 파싱합니다.
    """
    if not raw:
        return None, ""
    match = re.match(r"^([^(]+)\(([0-9.]+)\)$", raw)
    if match:
        pred = match.group(1).strip()
        try:
            score = float(match.group(2))
        except ValueError:
            score = None
        return score, pred
    # 점수 없이 예측만 있는 경우
    return None, raw


def extract_vep_annotations_from_vcf(
    vcf_path: str,
) -> Dict[str, Dict[str, Any]]:
    """
    VEP-annotated VCF 파일 전체를 읽어 각 변이의 CSQ annotation을 추출합니다.

    반환 딕셔너리 키: "CHROM:POS:REF:ALT" (예: "chr1:12345:A:G")
    반환 딕셔너리 값: _normalize_csq_record() 결과

    이 함수는 vcf_parser.py의 parse_vcf()에서 VEP VCF를 전처리할 때 사용됩니다.
    pysam이 없는 환경에서는 직접 파싱하는 fallback을 사용합니다.
    """
    annotations: Dict[str, Dict[str, Any]] = {}
    csq_fields: Optional[List[str]] = None

    try:
        import pysam  # type: ignore

        with pysam.VariantFile(vcf_path) as vcf:
            # CSQ 필드 목록을 헤더에서 추출
            if "CSQ" in vcf.header.info:
                csq_desc = vcf.header.info["CSQ"].description or ""
                idxfmt = csq_desc.lower().find("format:")
                if idxfmt >= 0:
                    tail = csq_desc[idxfmt + len("format:") :].strip()
                    if tail.startswith('"'):
                        tail = tail[1:]
                    fields_part = tail.split('"')[0].strip()
                    if fields_part:
                        csq_fields = [f.strip() for f in fields_part.split("|")]

            if not csq_fields:
                logger.warning(
                    "VEP CSQ header not found in %s; using default field list", vcf_path
                )
                csq_fields = _DEFAULT_CSQ_FIELDS

            for rec in vcf.fetch():
                csq_raw = rec.info.get("CSQ")
                if not csq_raw:
                    continue
                # pysam은 CSQ를 튜플로 반환
                csq_str = ",".join(csq_raw) if isinstance(csq_raw, tuple) else str(csq_raw)
                for alt in rec.alts or []:
                    key = f"{rec.chrom}:{rec.pos}:{rec.ref}:{alt}"
                    parsed = parse_csq_record(csq_str, csq_fields, alt_allele=alt)
                    if parsed:
                        annotations[key] = parsed

    except ImportError:
        logger.warning("pysam not available; falling back to line-by-line VCF parsing")
        annotations = _extract_vep_annotations_fallback(vcf_path)
    except Exception as e:
        logger.error("Failed to extract VEP annotations from %s: %s", vcf_path, e)

    logger.info(
        "VEP CSQ extraction complete: %d variants annotated from %s",
        len(annotations),
        vcf_path,
    )
    return annotations


def _extract_vep_annotations_fallback(vcf_path: str) -> Dict[str, Dict[str, Any]]:
    """
    pysam 없이 VCF를 직접 파싱하는 fallback 구현.
    bgzip 압축 파일(.vcf.gz)은 gzip으로 읽습니다.
    """
    import gzip

    annotations: Dict[str, Dict[str, Any]] = {}
    csq_fields: Optional[List[str]] = None
    header_lines: List[str] = []

    opener = gzip.open if vcf_path.endswith(".gz") else open

    try:
        with opener(vcf_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    header_lines.append(line)
                    continue
                if line.startswith("#CHROM"):
                    # 헤더 파싱
                    csq_fields = parse_csq_header("\n".join(header_lines))
                    if not csq_fields:
                        csq_fields = _DEFAULT_CSQ_FIELDS
                    continue

                parts = line.split("\t")
                if len(parts) < 8:
                    continue

                chrom, pos, _, ref, alt_field, _, _, info_field = parts[:8]
                alts = alt_field.split(",")

                # INFO 필드에서 CSQ 추출
                csq_str = ""
                for token in info_field.split(";"):
                    if token.startswith("CSQ="):
                        csq_str = token[4:]
                        break

                if not csq_str:
                    continue

                for alt in alts:
                    key = f"{chrom}:{pos}:{ref}:{alt}"
                    parsed = parse_csq_record(csq_str, csq_fields or _DEFAULT_CSQ_FIELDS, alt_allele=alt)
                    if parsed:
                        annotations[key] = parsed

    except Exception as e:
        logger.error("Fallback VEP parsing failed for %s: %s", vcf_path, e)

    return annotations


def is_vep_annotated_vcf(vcf_path: str) -> bool:
    """
    VCF 파일이 VEP annotation을 포함하는지 헤더를 확인합니다.
    service-daemon의 plugin.py에서 VEP vs snpEff 분기 판단에 사용합니다.
    """
    import gzip

    opener = gzip.open if vcf_path.endswith(".gz") else open
    try:
        with opener(vcf_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                # Any CSQ INFO (VEP, Funcotator, or other) — skip redundant snpEff and try CSQ pre-parse
                if "##INFO=<ID=CSQ" in line:
                    return True
    except Exception:
        pass
    return False
