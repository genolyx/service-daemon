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

    # ── Priority 1: Frequency Rules (BA1 / BS1) – immediate return (원본 정렬) ──
    if gnomad_af is not None:
        try:
            af_val = float(gnomad_af)
            if af_val > 0.05:
                return {
                    "classification": "Benign",
                    "criteria_met": ["BA1"],
                    "reasoning": f"Common variant (gnomAD AF={af_val:.4f})",
                    "confidence": "high",
                }
            if af_val > 0.01:
                return {
                    "classification": "Likely Benign",
                    "criteria_met": ["BS1"],
                    "reasoning": f"Elevated allele frequency (gnomAD AF={af_val:.4f})",
                    "confidence": "high",
                }
        except (TypeError, ValueError):
            pass

    # ── Priority 2: ClinVar Rules – immediate return (원본 정렬) ─────────────
    if clinvar_sig:
        if "pathogenic" in clinvar_sig:
            # "Pathogenic/Likely pathogenic" 는 Pathogenic 으로
            if "likely" in clinvar_sig and "pathogenic/likely" not in clinvar_sig:
                return {
                    "classification": "Likely Pathogenic",
                    "criteria_met": ["PP5"],
                    "reasoning": f"ClinVar: {clinvar_sig}",
                    "confidence": "high" if clinvar_stars >= 2 else "medium",
                }
            return {
                "classification": "Pathogenic",
                "criteria_met": ["PS1", "PP5"],
                "reasoning": f"ClinVar Pathogenic ({clinvar_stars} stars)",
                "confidence": "high" if clinvar_stars >= 2 else "medium",
            }
        if "benign" in clinvar_sig:
            if "likely" in clinvar_sig and "benign/likely" not in clinvar_sig:
                return {
                    "classification": "Likely Benign",
                    "criteria_met": ["BP6"],
                    "reasoning": f"ClinVar: {clinvar_sig}",
                    "confidence": "high" if clinvar_stars >= 2 else "medium",
                }
            return {
                "classification": "Benign",
                "criteria_met": ["BP6", "BA1" if gnomad_af and gnomad_af > 0.05 else "BP6"],
                "reasoning": f"ClinVar Benign ({clinvar_stars} stars)",
                "confidence": "high" if clinvar_stars >= 2 else "medium",
            }

    # ── Priority 3: Variant Effect Rules ─────────────────────────────────────

    # PVS1: LOF variant
    if effects_set & LOF_EFFECTS:
        criteria.append("PVS1")
        reasoning_parts.append(f"Loss-of-function variant ({effect})")

    # PM2: Absent / rare in population
    if gnomad_af is not None:
        try:
            af_val = float(gnomad_af)
            if af_val < 0.0001:
                criteria.append("PM2")
                reasoning_parts.append(f"Rare in population (gnomAD AF={af_val:.6f})")
            elif af_val > 0.001:
                criteria.append("BS2")
                reasoning_parts.append(f"Observed in healthy controls (gnomAD AF={af_val:.4f})")
        except (TypeError, ValueError):
            pass
    else:
        criteria.append("PM2")
        reasoning_parts.append("Not found in gnomAD")

    # PM2_Supporting: moderate functional impact
    if effects_set & MODERATE_EFFECTS:
        criteria.append("PM2_Supporting")
        reasoning_parts.append(f"Moderate functional impact ({effect})")

    # BP7: silent / non-coding only
    if (effects_set & BENIGN_EFFECTS
            and not (effects_set & LOF_EFFECTS)
            and not (effects_set & MODERATE_EFFECTS)):
        criteria.append("BP7")
        reasoning_parts.append(f"Silent/non-coding variant ({effect})")

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

    # Stand-alone Benign (BA1 already handled in classify_acmg_lite early-return,
    # but keep here as safety net for direct calls to _determine_classification)
    if ba:
        return "Benign"

    # Benign (two strong benign, or strong + supporting)
    if len(bs) >= 2:
        return "Benign"
    if bs and bp:
        return "Likely Benign"
    if len(bp) >= 2:
        return "Likely Benign"

    # Pathogenic
    if pvs and (ps or len(pm) >= 1):
        return "Pathogenic"
    if len(ps) >= 2:
        return "Pathogenic"
    if ps and (len(pm) >= 1 or len(pp) >= 2):
        return "Pathogenic"

    # Likely Pathogenic
    if pvs and (pm or pp):
        return "Likely Pathogenic"
    if ps and (pm or pp):
        return "Likely Pathogenic"
    if len(pm) >= 3:
        return "Likely Pathogenic"
    if len(pm) >= 2 and len(pp) >= 2:
        return "Likely Pathogenic"

    # VUS (default when evidence is mixed or insufficient)
    if pvs or ps or pm or pp:
        return "VUS"

    # Likely benign (single benign evidence)
    if bs or bp:
        return "Likely Benign"

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
    "classification": "Pathogenic|Likely Pathogenic|VUS|Likely Benign|Benign",
    "criteria_met": ["PVS1", "PM2", ...],
    "reasoning": "Detailed reasoning...",
    "confidence": "high|medium|low"
}
"""


async def classify_acmg_ai(
    variant: Dict[str, Any],
    api_key: str = "",
    model: str = "",
    provider: str = "",
) -> Optional[Dict[str, Any]]:
    """
    AI 기반 ACMG 분류 (Gemini / OpenAI 지원).

    provider/model은 settings에서 읽어 사용하므로 직접 전달 불필요.
    """
    import os
    from .._resolve_ai import call_ai_json

    if not provider:
        try:
            from ...config import settings
            provider = settings.acmg_ai_provider
        except Exception:
            provider = os.environ.get("ACMG_AI_PROVIDER", "gemini")

    if not model:
        try:
            from ...config import settings
            model = settings.acmg_ai_model
        except Exception:
            model = "gemini-2.0-flash"

    if not api_key:
        try:
            from ...config import settings
            api_key = settings.gemini_api_key if provider == "gemini" else settings.acmg_ai_api_key
        except Exception:
            pass
    if not api_key:
        api_key = os.environ.get("GEMINI_API_KEY" if provider == "gemini" else "OPENAI_API_KEY", "")

    if not api_key:
        logger.warning(f"AI API key not configured for provider={provider}, skipping AI classification")
        return None

    variant_text = _format_variant_for_ai(variant)
    result = await call_ai_json(
        provider=provider,
        model=model,
        api_key=api_key,
        system_prompt=ACMG_AI_SYSTEM_PROMPT,
        user_content=variant_text,
    )
    if not result or "classification" not in result:
        return None
    return {
        "classification": result.get("classification", "VUS"),
        "criteria_met": result.get("criteria_met", []),
        "reasoning": result.get("reasoning", ""),
        "confidence": result.get("confidence", "low"),
        "source": f"AI({provider})",
    }


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
    provider: str = "",
    literature: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Rule-based + (선택적) Literature-enhanced AI 결합 분류.

    흐름:
      1. Rule-based (항상 실행)
      2. 결과가 VUS이고 literature 있으면 → Literature-enhanced AI로 강화 시도
      3. VUS가 아닌 경우(B/LB/LP/P)에는 AI 호출 불필요 (이미 명확)

    Args:
        variant:     annotator.annotate() 결과 딕셔너리
        use_ai:      AI 분류 활성화 여부
        api_key:     AI API 키 (Gemini 또는 OpenAI)
        provider:    "gemini" | "openai" (기본값: settings에서 읽음)
        literature:  search_variant_literature() 결과 (없으면 AI는 문헌 없이 동작)

    Returns:
        {
            "rule_based": {...},
            "ai": {...} or None,
            "literature_enhanced": bool,
            "final_classification": str,
            "final_criteria": list,
            "final_reasoning": str,
        }
    """
    rule_result = classify_acmg_lite(variant)

    ai_result = None
    literature_enhanced = False

    if use_ai and rule_result["classification"] == "VUS":
        if literature and literature.get("articles"):
            ai_result = await classify_acmg_with_literature(
                variant, literature["articles"],
                api_key=api_key, provider=provider,
            )
            literature_enhanced = True
        else:
            ai_result = await classify_acmg_ai(
                variant, api_key=api_key, provider=provider,
            )

    # 최종 분류: AI confidence=high이면 AI 우선, 아니면 rule-based
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
        "literature_enhanced": literature_enhanced,
        "final_classification": final_class,
        "final_criteria": final_criteria,
        "final_reasoning": final_reasoning,
    }


# ══════════════════════════════════════════════════════════════
# Literature-enhanced VUS Classification
# ══════════════════════════════════════════════════════════════

LITERATURE_ACMG_SYSTEM_PROMPT = """You are a clinical genetics expert specializing in ACMG/AMP variant classification.
You are given a variant with its annotation data AND relevant PubMed abstracts.

Use the literature to identify additional ACMG evidence criteria:
  - PS3/BS3: Functional studies showing damaging or benign effect
  - PP4:     Phenotype matches disease known to be caused by this gene
  - PM1:     Located in mutational hotspot or critical functional domain
  - PP3/BP4: Computational predictions (mention if cited in literature)

IMPORTANT:
  - Only cite criteria that are DIRECTLY supported by the provided abstracts
  - Do NOT invent or hallucinate criteria
  - If abstracts are irrelevant, return the same classification as rule_based
  - Cite PMIDs when applying PS3/PP4

Respond ONLY in JSON:
{
    "classification": "Pathogenic|Likely Pathogenic|VUS|Likely Benign|Benign",
    "criteria_met": ["PS3:PMID_12345678", "PP4", ...],
    "reasoning": "Brief reasoning with PMID citations",
    "confidence": "high|medium|low"
}"""


async def classify_acmg_with_literature(
    variant: Dict[str, Any],
    articles: List[Dict[str, Any]],
    api_key: str = "",
    model: str = "",
    provider: str = "",
) -> Optional[Dict[str, Any]]:
    """
    PubMed 논문 초록을 컨텍스트로 제공하여 VUS를 추가 분류합니다.
    Rule-based에서 VUS가 나왔을 때만 호출됩니다.
    """
    import os
    from .._resolve_ai import call_ai_json

    if not provider:
        try:
            from ...config import settings
            provider = settings.acmg_ai_provider
        except Exception:
            provider = os.environ.get("ACMG_AI_PROVIDER", "gemini")

    if not model:
        try:
            from ...config import settings
            model = settings.acmg_ai_model
        except Exception:
            model = "gemini-2.0-flash"

    if not api_key:
        try:
            from ...config import settings
            api_key = settings.gemini_api_key if provider == "gemini" else settings.acmg_ai_api_key
        except Exception:
            pass
    if not api_key:
        api_key = os.environ.get("GEMINI_API_KEY" if provider == "gemini" else "OPENAI_API_KEY", "")

    if not api_key:
        logger.warning(f"AI API key not configured for provider={provider}, skipping literature classification")
        return None

    top_articles = sorted(articles, key=lambda a: a.get("relevance_score", 0), reverse=True)[:5]
    if not top_articles:
        return None

    abstracts_text = "\n\n".join(
        f"PMID {a.get('pmid', '?')} ({a.get('pub_date', '')}): {a.get('title', '')}\n"
        f"{a.get('abstract', '')[:800]}"
        for a in top_articles
    )
    variant_text = _format_variant_for_ai(variant)
    user_content = (
        f"VARIANT:\n{variant_text}\n\n"
        f"RULE-BASED RESULT: VUS (no definitive criteria met)\n\n"
        f"RELEVANT LITERATURE ({len(top_articles)} papers):\n{abstracts_text}\n\n"
        "Based on this literature, can you identify any additional ACMG criteria? "
        "If the literature is not relevant to this specific variant, keep classification as VUS."
    )

    result = await call_ai_json(
        provider=provider,
        model=model,
        api_key=api_key,
        system_prompt=LITERATURE_ACMG_SYSTEM_PROMPT,
        user_content=user_content,
    )
    if not result or "classification" not in result:
        return None
    return {
        "classification": result.get("classification", "VUS"),
        "criteria_met": result.get("criteria_met", []),
        "reasoning": result.get("reasoning", ""),
        "confidence": result.get("confidence", "low"),
        "source": f"AI+Literature({provider})",
        "pmids_used": [a.get("pmid") for a in top_articles if a.get("pmid")],
    }
