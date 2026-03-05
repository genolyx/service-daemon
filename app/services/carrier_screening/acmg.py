"""
ACMG Variant Classification Module

ACMG/AMP 가이드라인에 기반한 변이 분류를 수행합니다.
genetic_reporter/acmg_classifier.py의 로직을 모듈화한 것입니다.

두 가지 모드:
    1. Rule-based (Lite): ClinVar, gnomAD, 기능적 영향 기반 자동 분류
    2. AI-assisted: OpenAI API를 통한 심층 분류 (선택적)
"""

import re
import json
import logging
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# ACMG Rule-based Classification (Lite)
# ══════════════════════════════════════════════════════════════

# 기능적 영향 분류
LOF_EFFECTS = {
    "frameshift_variant", "stop_gained", "splice_acceptor_variant",
    "splice_donor_variant", "start_lost", "stop_lost",
}

MODERATE_EFFECTS = {
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "splice_region_variant",
}

BENIGN_EFFECTS = {
    "synonymous_variant", "intron_variant", "intergenic_variant",
    "upstream_gene_variant", "downstream_gene_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
}


def classify_acmg_lite(variant: Dict[str, Any]) -> Dict[str, Any]:
    """
    Rule-based ACMG 분류를 수행합니다.

    Args:
        variant: annotator.annotate()의 결과 딕셔너리

    Returns:
        {
            "classification": str,  # Pathogenic / Likely pathogenic / VUS / Likely benign / Benign
            "criteria_met": list,   # 적용된 ACMG 기준 코드 목록
            "reasoning": str,       # 분류 근거 요약
            "confidence": str,      # high / medium / low
        }
    """
    criteria = []
    reasoning_parts = []

    effect = (variant.get("effect") or "").strip()
    effects_set = set(re.split(r"[&,]", effect)) if effect else set()

    clinvar_sig = (variant.get("clinvar_sig_primary") or "").strip().lower()
    clinvar_stars = int(variant.get("clinvar_stars") or 0)
    gnomad_af = variant.get("gnomad_af")
    zygosity = (variant.get("zygosity") or "").strip()

    # ── Pathogenic Evidence ──────────────────────────────

    # PVS1: Null variant (LOF) in a gene where LOF is a known mechanism
    if effects_set & LOF_EFFECTS:
        criteria.append("PVS1")
        reasoning_parts.append(f"Loss-of-function variant ({effect})")

    # PS1: Same amino acid change as established pathogenic variant
    if clinvar_sig in ("pathogenic",) and clinvar_stars >= 2:
        criteria.append("PS1")
        reasoning_parts.append(f"ClinVar Pathogenic ({clinvar_stars} stars)")

    # PM1: Located in a mutational hot spot / functional domain
    # (simplified: moderate effect in known gene)
    if effects_set & MODERATE_EFFECTS:
        criteria.append("PM2_Supporting")
        reasoning_parts.append(f"Moderate functional impact ({effect})")

    # PM2: Absent from controls (gnomAD AF < 0.0001)
    if gnomad_af is not None and gnomad_af < 0.0001:
        criteria.append("PM2")
        reasoning_parts.append(f"Rare in population (gnomAD AF={gnomad_af:.6f})")
    elif gnomad_af is None:
        criteria.append("PM2")
        reasoning_parts.append("Not found in gnomAD")

    # PP3: Computational evidence (missense in LOF-intolerant gene)
    # Simplified: ClinVar Likely pathogenic
    if clinvar_sig in ("likely pathogenic",):
        criteria.append("PP3")
        reasoning_parts.append("ClinVar Likely pathogenic")

    # PP5: Reputable source reports as pathogenic
    if clinvar_sig in ("pathogenic",) and clinvar_stars >= 1:
        criteria.append("PP5")
        reasoning_parts.append("ClinVar Pathogenic report")

    # ── Benign Evidence ──────────────────────────────────

    # BA1: Allele frequency > 5%
    if gnomad_af is not None and gnomad_af > 0.05:
        criteria.append("BA1")
        reasoning_parts.append(f"Common variant (gnomAD AF={gnomad_af:.4f})")

    # BS1: Allele frequency > expected for disorder (> 1%)
    elif gnomad_af is not None and gnomad_af > 0.01:
        criteria.append("BS1")
        reasoning_parts.append(f"Elevated frequency (gnomAD AF={gnomad_af:.4f})")

    # BS2: Observed in healthy adult (gnomAD AF > 0.001)
    elif gnomad_af is not None and gnomad_af > 0.001:
        criteria.append("BS2")
        reasoning_parts.append(f"Observed in healthy controls (gnomAD AF={gnomad_af:.4f})")

    # BP1: Missense in gene where only truncating cause disease
    if effects_set & BENIGN_EFFECTS and not (effects_set & LOF_EFFECTS) and not (effects_set & MODERATE_EFFECTS):
        criteria.append("BP7")
        reasoning_parts.append(f"Silent/non-coding variant ({effect})")

    # BP6: Reputable source reports as benign
    if clinvar_sig in ("benign", "likely benign"):
        criteria.append("BP6")
        reasoning_parts.append(f"ClinVar {clinvar_sig}")

    # ── Classification Logic ─────────────────────────────

    classification = _determine_classification(criteria)
    confidence = _determine_confidence(criteria, clinvar_stars)

    return {
        "classification": classification,
        "criteria_met": criteria,
        "reasoning": "; ".join(reasoning_parts) if reasoning_parts else "No specific criteria met",
        "confidence": confidence,
    }


def _determine_classification(criteria: List[str]) -> str:
    """ACMG 기준 조합에 따른 최종 분류를 결정합니다."""
    has = set(criteria)

    # Pathogenic: PVS1 + PM2, or PS1 + PM2, etc.
    pvs = [c for c in has if c.startswith("PVS")]
    ps = [c for c in has if c.startswith("PS")]
    pm = [c for c in has if c.startswith("PM")]
    pp = [c for c in has if c.startswith("PP")]
    ba = [c for c in has if c.startswith("BA")]
    bs = [c for c in has if c.startswith("BS")]
    bp = [c for c in has if c.startswith("BP")]

    # Stand-alone Benign
    if ba:
        return "Benign"

    # Benign
    if len(bs) >= 2:
        return "Benign"
    if bs and bp:
        return "Likely benign"
    if len(bp) >= 2:
        return "Likely benign"

    # Pathogenic
    if pvs and (ps or len(pm) >= 1):
        return "Pathogenic"
    if len(ps) >= 2:
        return "Pathogenic"
    if ps and (len(pm) >= 1 or len(pp) >= 2):
        return "Pathogenic"

    # Likely pathogenic
    if pvs and pm:
        return "Likely pathogenic"
    if pvs and pp:
        return "Likely pathogenic"
    if ps and pm:
        return "Likely pathogenic"
    if len(pm) >= 3:
        return "Likely pathogenic"
    if len(pm) >= 2 and len(pp) >= 2:
        return "Likely pathogenic"

    # VUS (default when evidence is mixed or insufficient)
    if pvs or ps or pm or pp:
        return "VUS"

    # Likely benign (single benign evidence)
    if bs or bp:
        return "Likely benign"

    return "VUS"


def _determine_confidence(criteria: List[str], clinvar_stars: int) -> str:
    """분류 신뢰도를 결정합니다."""
    if clinvar_stars >= 3:
        return "high"
    if clinvar_stars >= 2 and len(criteria) >= 2:
        return "high"
    if len(criteria) >= 3:
        return "medium"
    if len(criteria) >= 1:
        return "low"
    return "low"


# ══════════════════════════════════════════════════════════════
# ACMG AI-assisted Classification (Optional)
# ══════════════════════════════════════════════════════════════

ACMG_AI_SYSTEM_PROMPT = """You are a clinical genetics expert specializing in ACMG/AMP variant classification.
Given a variant with its annotation data, classify it according to ACMG/AMP guidelines.

Respond in JSON format:
{
    "classification": "Pathogenic|Likely pathogenic|VUS|Likely benign|Benign",
    "criteria_met": ["PVS1", "PM2", ...],
    "reasoning": "Detailed reasoning...",
    "confidence": "high|medium|low"
}
"""


async def classify_acmg_ai(
    variant: Dict[str, Any],
    api_key: str = "",
    model: str = "gpt-4.1-mini",
) -> Optional[Dict[str, Any]]:
    """
    OpenAI API를 사용한 AI 기반 ACMG 분류.

    Args:
        variant: annotator.annotate()의 결과 딕셔너리
        api_key: OpenAI API 키
        model: 사용할 모델

    Returns:
        분류 결과 딕셔너리 또는 None (실패 시)
    """
    if not api_key:
        import os
        api_key = os.environ.get("OPENAI_API_KEY", "")
    if not api_key:
        logger.warning("OpenAI API key not configured, skipping AI classification")
        return None

    try:
        from openai import AsyncOpenAI

        # 변이 정보를 프롬프트용 텍스트로 변환
        variant_text = _format_variant_for_ai(variant)

        client = AsyncOpenAI(api_key=api_key)
        response = await client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": ACMG_AI_SYSTEM_PROMPT},
                {"role": "user", "content": variant_text},
            ],
            temperature=0.1,
            max_tokens=1000,
            response_format={"type": "json_object"},
        )

        content = response.choices[0].message.content
        result = json.loads(content)

        # 필수 필드 검증
        if "classification" not in result:
            logger.warning("AI classification missing 'classification' field")
            return None

        return {
            "classification": result.get("classification", "VUS"),
            "criteria_met": result.get("criteria_met", []),
            "reasoning": result.get("reasoning", ""),
            "confidence": result.get("confidence", "low"),
            "source": "AI",
        }

    except Exception as e:
        logger.error(f"AI ACMG classification failed: {e}")
        return None


def _format_variant_for_ai(variant: Dict[str, Any]) -> str:
    """변이 정보를 AI 프롬프트용 텍스트로 변환합니다."""
    lines = [
        f"Gene: {variant.get('gene', 'Unknown')}",
        f"Position: {variant.get('chrom', '')}:{variant.get('pos', '')}",
        f"Change: {variant.get('ref', '')} > {variant.get('alt', '')}",
        f"HGVS.c: {variant.get('hgvsc', '')}",
        f"HGVS.p: {variant.get('hgvsp', '')}",
        f"Effect: {variant.get('effect', '')}",
        f"Zygosity: {variant.get('zygosity', '')}",
        f"gnomAD AF: {variant.get('gnomad_af', 'Not found')}",
        f"ClinVar: {variant.get('clinvar_sig_primary', 'Not found')} "
        f"(stars: {variant.get('clinvar_stars', 0)})",
        f"ClinVar disease: {variant.get('clinvar_dn', '')}",
        f"ClinGen HI score: {variant.get('clingen_hi_score', 'N/A')}",
        f"ClinGen TS score: {variant.get('clingen_ts_score', 'N/A')}",
        f"Depth: {variant.get('dp', 'N/A')}",
        f"VAF: {variant.get('vaf', 'N/A')}",
    ]
    return "\n".join(lines)


# ══════════════════════════════════════════════════════════════
# Combined Classification
# ══════════════════════════════════════════════════════════════

async def classify_variant(
    variant: Dict[str, Any],
    use_ai: bool = False,
    api_key: str = "",
) -> Dict[str, Any]:
    """
    Rule-based + AI (선택적) 결합 분류.

    Returns:
        {
            "rule_based": {...},
            "ai": {...} or None,
            "final_classification": str,
            "final_criteria": list,
            "final_reasoning": str,
        }
    """
    # Rule-based 분류
    rule_result = classify_acmg_lite(variant)

    # AI 분류 (선택적)
    ai_result = None
    if use_ai:
        ai_result = await classify_acmg_ai(variant, api_key=api_key)

    # 최종 분류 결정
    if ai_result and ai_result.get("confidence") == "high":
        final_class = ai_result["classification"]
        final_criteria = ai_result.get("criteria_met", [])
        final_reasoning = ai_result.get("reasoning", "")
    else:
        final_class = rule_result["classification"]
        final_criteria = rule_result["criteria_met"]
        final_reasoning = rule_result["reasoning"]

    return {
        "rule_based": rule_result,
        "ai": ai_result,
        "final_classification": final_class,
        "final_criteria": final_criteria,
        "final_reasoning": final_reasoning,
    }
