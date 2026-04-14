"""
PGx (PharmCAT) pipeline outputs under <analysis_or_output>/pgx/ → result.json + customer PDF.

Expected layout (gx-exome):
  pgx/pgx_meta.json
  pgx/pgx_summary.txt
  pgx/pgx_result.json          (artifact pointers; may be large — not embedded wholesale)
  pgx/<sample>_pgx.report.html (full PharmCAT HTML — not embedded in PDF; use summary + meta)
  pgx/<sample>_pgx.report.json (PharmCAT reporter JSON — drug recommendations)

The customer PDF uses the plain-text summary (human-readable, WeasyPrint-safe).

Portal also gets ``gene_results``: one row per gene from ``pgx_result.json`` → ``phenotype``
(PharmCAT ``geneReports``), with ``reviewer_confirmed`` / ``reviewer_comment`` merged on reprocess.

Drug recommendations are extracted from the PharmCAT reporter JSON (``*.report.json``)
or from the ``pgx_result.json`` bundle (if it contains a ``report`` or ``drugs`` key).
Only actionable gene/drug pairs are included on the customer PDF.
"""

from __future__ import annotations

import csv
import glob
import html
import io
import json
import logging
import os
import re
from collections import Counter
from typing import AbstractSet, Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Mirrors gx-exome `pgx.nf` PYSUMMARY (actionable vs normal)
_RISK_FUNCTIONS = frozenset({"no function", "decreased function", "unfavorable response allele"})
_SKIP_PHENOTYPES = frozenset({"no result", "n/a", ""})


def extract_gene_results_from_phenotype(phenotype: Any) -> List[Dict[str, Any]]:
    """
    Build one row per gene from PharmCAT phenotype JSON (``geneReports`` / recommendation diplotypes).

    Handles both v2 format (nested ``geneReports.CPIC.<gene>`` / ``geneReports.DPWG.<gene>``)
    and v3 format (flat ``geneReports.<gene>``).
    """
    if not isinstance(phenotype, dict):
        return []
    gene_reports = phenotype.get("geneReports")
    if not isinstance(gene_reports, dict):
        return []

    # Detect format: v2 has CPIC/DPWG keys mapping to dicts of genes;
    # v3 has gene names mapping directly to gene data dicts.
    source_gene_pairs: List[tuple] = []
    if "CPIC" in gene_reports or "DPWG" in gene_reports:
        # v2 format
        for source in ("CPIC", "DPWG"):
            reports = gene_reports.get(source)
            if not isinstance(reports, dict):
                continue
            for gene_name in sorted(reports.keys()):
                source_gene_pairs.append((source, gene_name, reports[gene_name]))
    else:
        # v3 format — flat dict of gene_name -> gene_data
        for gene_name in sorted(gene_reports.keys()):
            gd = gene_reports[gene_name]
            if isinstance(gd, dict) and gd.get("recommendationDiplotypes"):
                src = (gd.get("alleleDefinitionSource") or "CPIC").upper()
                if src == "CLINPGX":
                    src = "CPIC"
                source_gene_pairs.append((src, gene_name, gd))

    rows: List[Dict[str, Any]] = []
    seen_genes: set = set()
    for source, gene_name, gene_data in source_gene_pairs:
        if gene_name in seen_genes:
            continue
        if not isinstance(gene_data, dict):
            continue
        call_src = (gene_data.get("callSource") or "").strip()
        rec_dips = gene_data.get("recommendationDiplotypes")
        if not isinstance(rec_dips, list):
            continue
        for dip in rec_dips:
            if not isinstance(dip, dict):
                continue
            a1 = dip.get("allele1") if isinstance(dip.get("allele1"), dict) else {}
            a2 = dip.get("allele2") if isinstance(dip.get("allele2"), dict) else {}
            n1 = (a1.get("name") or "").strip()
            n2 = (a2.get("name") or "").strip()
            fn1 = (a1.get("function") or "").strip()
            fn2 = (a2.get("function") or "").strip()
            phenotypes = dip.get("phenotypes")
            if not isinstance(phenotypes, list):
                phenotypes = []
            activity = dip.get("activityScore")
            diplotype = f"{n1}/{n2}" if n2 else n1
            phenotype_str = ", ".join(p for p in phenotypes if p) if phenotypes else ""
            low = phenotype_str.lower()
            if low in _SKIP_PHENOTYPES:
                continue
            functions = [f for f in (fn1, fn2) if f]
            has_risk = any((f or "").lower() in _RISK_FUNCTIONS for f in functions)
            cat = "actionable" if has_risk else "normal"
            rows.append(
                {
                    "gene": gene_name,
                    "guideline_source": source,
                    "diplotype": diplotype,
                    "phenotype": phenotype_str,
                    "activity_score": activity,
                    "allele1_function": fn1,
                    "allele2_function": fn2,
                    "call_source": call_src,
                    "category": cat,
                }
            )
            seen_genes.add(gene_name)
            break
    rows.sort(
        key=lambda r: (
            0 if r.get("category") == "actionable" else 1,
            str(r.get("gene") or ""),
        )
    )
    return rows


def extract_drug_recommendations(
    report_data: Any,
    actionable_genes: Optional[AbstractSet[str]] = None,
) -> List[Dict[str, Any]]:
    """
    Extract drug prescribing recommendations from PharmCAT reporter JSON.

    Handles multiple PharmCAT JSON layouts:
      - v2/v3 reporter JSON: top-level ``drugs`` array
      - bundled ``pgx_result.json``: ``report.drugs`` or top-level ``drugs``
      - v2 phenotype-embedded: ``geneReports.CPIC/DPWG.<gene>.recommendations``

    Returns a list of dicts with: drug, gene, source, classification, phenotype,
    implication, recommendation.  Sorted actionable-first, then by drug name.
    """
    if not isinstance(report_data, dict):
        return []

    drugs_array: Optional[List] = None

    if isinstance(report_data.get("drugs"), list):
        drugs_array = report_data["drugs"]
    elif isinstance(report_data.get("report"), dict):
        rpt = report_data["report"]
        if isinstance(rpt.get("drugs"), list):
            drugs_array = rpt["drugs"]
    if isinstance(report_data.get("groups"), list):
        drugs_array = drugs_array or []
        for grp in report_data["groups"]:
            if isinstance(grp, dict) and isinstance(grp.get("drugs"), list):
                for d in grp["drugs"]:
                    drugs_array.append(d)

    if not drugs_array:
        phen = report_data.get("phenotype")
        if isinstance(phen, dict):
            drugs_array = _drug_recs_from_phenotype_gene_reports(phen)

    if not drugs_array:
        return []

    actionable_up = (
        {g.upper() for g in actionable_genes} if actionable_genes else None
    )
    rows: List[Dict[str, Any]] = []

    for drug_obj in drugs_array:
        if not isinstance(drug_obj, dict):
            continue
        drug_name = ""
        if isinstance(drug_obj.get("name"), str):
            drug_name = drug_obj["name"].strip()
        elif isinstance(drug_obj.get("drug"), dict):
            drug_name = (drug_obj["drug"].get("name") or "").strip()
        elif isinstance(drug_obj.get("drug"), str):
            drug_name = drug_obj["drug"].strip()
        if not drug_name:
            continue

        annotations = drug_obj.get("guidelineAnnotations") or drug_obj.get("guidelines") or []
        if not isinstance(annotations, list):
            annotations = [annotations] if isinstance(annotations, dict) else []

        for ann in annotations:
            if not isinstance(ann, dict):
                continue
            source = (ann.get("source") or "").strip()
            classification = (ann.get("classification") or "").strip()
            genotype = (ann.get("genotype") or "").strip()

            gene = ""
            if genotype and ":" in genotype:
                gene = genotype.split(":")[0].strip()
            elif isinstance(ann.get("gene"), str):
                gene = ann["gene"].strip()
            elif isinstance(drug_obj.get("gene"), str):
                gene = drug_obj["gene"].strip()

            if actionable_up and gene and gene.upper() not in actionable_up:
                continue

            phenotype = ""
            if isinstance(ann.get("phenotype"), str):
                phenotype = ann["phenotype"].strip()
            elif isinstance(ann.get("phenotypes"), dict):
                phenotype = ", ".join(
                    f"{k}: {v}" for k, v in ann["phenotypes"].items() if v
                )

            implications_raw = ann.get("implications") or ann.get("implication") or ""
            if isinstance(implications_raw, dict):
                implication = "; ".join(str(v) for v in implications_raw.values() if v)
            elif isinstance(implications_raw, list):
                implication = "; ".join(str(v) for v in implications_raw if v)
            else:
                implication = str(implications_raw).strip()

            rec_raw = (
                ann.get("recommendation")
                or ann.get("drugRecommendation")
                or ann.get("drugRecommendations")
                or ""
            )
            if isinstance(rec_raw, list):
                recommendation = " ".join(
                    (r.get("text") if isinstance(r, dict) else str(r))
                    for r in rec_raw if r
                ).strip()
            elif isinstance(rec_raw, dict):
                recommendation = (rec_raw.get("text") or "").strip()
            else:
                recommendation = str(rec_raw).strip()

            if not recommendation and not implication:
                continue

            rows.append({
                "drug": drug_name,
                "gene": gene,
                "source": source,
                "classification": classification,
                "genotype": genotype,
                "phenotype": phenotype,
                "implication": implication,
                "recommendation": recommendation,
            })

    rows.sort(key=lambda r: (r.get("drug") or "", r.get("gene") or ""))
    return rows


def _drug_recs_from_phenotype_gene_reports(phenotype: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Fallback: build pseudo-drug objects from phenotype geneReports recommendations."""
    drugs: List[Dict[str, Any]] = []
    gene_reports = phenotype.get("geneReports")
    if not isinstance(gene_reports, dict):
        return drugs

    # Build (source, gene_name, gene_data) tuples — v2 vs v3 format detection
    source_gene_pairs: List[tuple] = []
    if "CPIC" in gene_reports or "DPWG" in gene_reports:
        for source in ("CPIC", "DPWG"):
            reports = gene_reports.get(source)
            if isinstance(reports, dict):
                for gn, gd in reports.items():
                    source_gene_pairs.append((source, gn, gd))
    else:
        for gn, gd in gene_reports.items():
            if isinstance(gd, dict):
                src = (gd.get("alleleDefinitionSource") or "CPIC").upper()
                if src == "CLINPGX":
                    src = "CPIC"
                source_gene_pairs.append((src, gn, gd))

    for source, gene_name, gene_data in source_gene_pairs:
        if not isinstance(gene_data, dict):
            continue
        recs = gene_data.get("recommendations")
        if not isinstance(recs, list):
            continue
        for rec in recs:
            if not isinstance(rec, dict):
                continue
            drug_info = rec.get("drug")
            dname = ""
            if isinstance(drug_info, dict):
                dname = (drug_info.get("name") or "").strip()
            elif isinstance(drug_info, str):
                dname = drug_info.strip()
            if not dname:
                continue
            drugs.append({
                "name": dname,
                "gene": gene_name,
                "guidelineAnnotations": [{
                    "source": source,
                    "gene": gene_name,
                    "classification": (rec.get("classification") or "").strip(),
                    "phenotype": rec.get("phenotype") or "",
                    "implications": rec.get("implications") or "",
                    "recommendation": rec.get("drugRecommendation") or rec.get("recommendation") or "",
                }],
            })
    return drugs


def _collect_drug_recommendations_from_pgx_dir(
    pgx_dir: str,
    actionable_genes: Optional[AbstractSet[str]] = None,
) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    """
    Scan ``pgx/`` for PharmCAT reporter JSON (``*.report.json``) and extract drug recs.

    Also checks ``pgx_result.json`` bundle for embedded report/drug data.
    Returns (drug_recommendations, source_file_basename).
    """
    recs: List[Dict[str, Any]] = []
    src: Optional[str] = None

    report_jsons = sorted(glob.glob(os.path.join(pgx_dir, "*.report.json")))
    for rj_path in report_jsons:
        try:
            with open(rj_path, "r", encoding="utf-8") as f:
                rj = json.load(f)
            extracted = extract_drug_recommendations(rj, actionable_genes)
            if extracted:
                recs = extracted
                src = os.path.basename(rj_path)
                break
        except Exception as e:
            logger.warning("[pgx] could not parse drug recs from %s: %s", rj_path, e)

    if not recs:
        bundle_path = os.path.join(pgx_dir, "pgx_result.json")
        if os.path.isfile(bundle_path):
            try:
                with open(bundle_path, "r", encoding="utf-8") as f:
                    bundle = json.load(f)
                extracted = extract_drug_recommendations(bundle, actionable_genes)
                if extracted:
                    recs = extracted
                    src = "pgx_result.json"
            except Exception as e:
                logger.warning("[pgx] could not parse drug recs from bundle: %s", e)

    return recs, src


_CPIC_DRUG_RECS: Dict[str, List[Dict[str, str]]] = {
    "CYP2C9": [
        {"drug": "Warfarin", "source": "CPIC", "implication": "Reduced warfarin metabolism; lower dose requirements.", "recommendation": "Consider reduced initial dose. Use validated dosing algorithms incorporating CYP2C9 genotype."},
        {"drug": "Celecoxib", "source": "CPIC", "implication": "Reduced celecoxib metabolism; increased exposure.", "recommendation": "Consider starting with 25–50% of the lowest recommended dose."},
        {"drug": "Phenytoin / Fosphenytoin", "source": "CPIC", "implication": "Reduced phenytoin clearance.", "recommendation": "Consider 25% reduction of starting dose; monitor drug levels closely."},
        {"drug": "NSAIDs (Flurbiprofen, Ibuprofen, Meloxicam, Piroxicam)", "source": "CPIC", "implication": "Reduced NSAID metabolism; increased risk of GI bleeding.", "recommendation": "Use lowest effective dose and shortest duration; consider alternative analgesics."},
    ],
    "CYP2D6": [
        {"drug": "Codeine", "source": "CPIC", "implication": "Reduced conversion of codeine to morphine; diminished analgesic effect.", "recommendation": "Consider alternative analgesic not metabolized by CYP2D6 (e.g. morphine, non-opioid)."},
        {"drug": "Tramadol", "source": "CPIC", "implication": "Reduced formation of active metabolite; diminished analgesic effect.", "recommendation": "Consider alternative analgesic; monitor for inadequate pain control."},
        {"drug": "Tamoxifen", "source": "CPIC", "implication": "Reduced formation of active metabolite endoxifen.", "recommendation": "Consider alternative endocrine therapy (e.g. aromatase inhibitor) or higher tamoxifen dose with monitoring."},
        {"drug": "Amitriptyline / Nortriptyline", "source": "CPIC", "implication": "Altered tricyclic antidepressant metabolism.", "recommendation": "Consider 25% dose reduction or alternative antidepressant; monitor drug levels."},
        {"drug": "Ondansetron / Tropisetron", "source": "CPIC", "implication": "Reduced conversion to active metabolite.", "recommendation": "Consider alternative antiemetic (e.g. granisetron)."},
    ],
    "CYP3A5": [
        {"drug": "Tacrolimus", "source": "CPIC", "implication": "CYP3A5 non-expresser; standard tacrolimus metabolism.", "recommendation": "Initiate at standard recommended dose (0.3 mg/kg/day). Adjust based on therapeutic drug monitoring."},
    ],
    "SLCO1B1": [
        {"drug": "Simvastatin", "source": "CPIC", "implication": "Increased simvastatin acid exposure; elevated myopathy risk.", "recommendation": "Prescribe ≤20 mg/day simvastatin or consider alternative statin (rosuvastatin, pravastatin, fluvastatin)."},
        {"drug": "Atorvastatin", "source": "CPIC", "implication": "Moderately increased atorvastatin exposure.", "recommendation": "Prescribe ≤40 mg/day; monitor for myopathy symptoms."},
        {"drug": "Rosuvastatin", "source": "CPIC", "implication": "Moderately increased rosuvastatin exposure.", "recommendation": "Prescribe ≤20 mg/day; monitor for muscle-related adverse effects."},
    ],
    "UGT1A1": [
        {"drug": "Irinotecan", "source": "CPIC", "implication": "Reduced UGT1A1 activity; increased SN-38 exposure and toxicity risk.", "recommendation": "Reduce starting dose by ≥ one level; monitor closely for neutropenia and diarrhea."},
        {"drug": "Atazanavir", "source": "CPIC", "implication": "Increased risk of hyperbilirubinemia (jaundice).", "recommendation": "Consider alternative antiretroviral if jaundice develops."},
    ],
    "CYP2C19": [
        {"drug": "Clopidogrel", "source": "CPIC", "implication": "Reduced formation of active metabolite; diminished antiplatelet effect.", "recommendation": "Consider alternative antiplatelet agent (e.g. prasugrel, ticagrelor)."},
        {"drug": "Voriconazole", "source": "CPIC", "implication": "Altered voriconazole exposure.", "recommendation": "Adjust dose per CPIC guidelines; consider therapeutic drug monitoring."},
        {"drug": "Escitalopram / Citalopram", "source": "CPIC", "implication": "Altered SSRI metabolism.", "recommendation": "Dose adjustment or alternative SSRI per CPIC guidelines."},
        {"drug": "Proton Pump Inhibitors (Omeprazole, Lansoprazole)", "source": "CPIC", "implication": "Altered PPI metabolism.", "recommendation": "Adjust dose per CPIC guidelines; consider therapeutic drug monitoring."},
    ],
    "DPYD": [
        {"drug": "Fluoropyrimidines (5-FU, Capecitabine)", "source": "CPIC", "implication": "Reduced DPD activity; increased risk of severe/fatal toxicity.", "recommendation": "Reduce dose by 50% for intermediate metabolizers; fluoropyrimidines are contraindicated for poor metabolizers."},
    ],
    "TPMT": [
        {"drug": "Thiopurines (Azathioprine, Mercaptopurine, Thioguanine)", "source": "CPIC", "implication": "Altered thiopurine metabolism; risk of myelosuppression.", "recommendation": "Reduce starting dose per CPIC guidelines; monitor blood counts closely."},
    ],
    "NUDT15": [
        {"drug": "Thiopurines (Azathioprine, Mercaptopurine, Thioguanine)", "source": "CPIC", "implication": "Increased risk of thiopurine-induced leukopenia.", "recommendation": "Reduce starting dose per CPIC guidelines; monitor blood counts closely."},
    ],
    "CYP2B6": [
        {"drug": "Efavirenz", "source": "CPIC", "implication": "Altered efavirenz metabolism.", "recommendation": "Dose adjustment per CPIC guidelines; consider alternative antiretroviral."},
    ],
    "VKORC1": [
        {"drug": "Warfarin", "source": "DPWG", "implication": "Altered warfarin sensitivity based on VKORC1 genotype.", "recommendation": "Use pharmacogenomics-based dosing algorithm incorporating VKORC1 genotype."},
    ],
}


def generate_cpic_drug_recommendations(
    gene_results: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Generate drug recommendations from built-in CPIC reference data based on gene results.

    Used as fallback when PharmCAT reporter JSON (with full drug annotations) is unavailable.
    Only generates recommendations for actionable genes.
    """
    recs: List[Dict[str, Any]] = []
    for row in gene_results:
        if not isinstance(row, dict):
            continue
        if row.get("category") != "actionable":
            continue
        gene = (row.get("gene") or "").strip()
        phenotype = (row.get("phenotype") or "").strip()
        if not gene:
            continue
        cpic_entries = _CPIC_DRUG_RECS.get(gene, [])
        for entry in cpic_entries:
            recs.append({
                "drug": entry["drug"],
                "gene": gene,
                "source": entry["source"],
                "classification": "CPIC Guideline",
                "phenotype": phenotype,
                "implication": entry["implication"],
                "recommendation": entry["recommendation"],
            })
    recs.sort(key=lambda r: (r.get("gene") or "", r.get("drug") or ""))
    return recs


def _read_pgx_custom_json(pgx_dir: str) -> Optional[Dict[str, Any]]:
    """Read ``pgx_custom_result.json`` if present (extended PGx / APOE tag SNPs)."""
    path = os.path.join(pgx_dir, "pgx_custom_result.json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        logger.warning("[pgx] could not read custom result: %s", e)
        return None
    return data if isinstance(data, dict) and not data.get("error") else None


def _custom_gene_rows_from_pgx_dict(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Flatten ``genes`` from pgx_custom_result.json into portal rows."""
    genes = data.get("genes")
    if not isinstance(genes, dict):
        return []
    rows: List[Dict[str, Any]] = []
    for gene_name in sorted(genes.keys()):
        for v in genes[gene_name]:
            if not isinstance(v, dict):
                continue
            zyg = (v.get("zygosity") or "").strip()
            is_variant = v.get("status") == "variant_found" and zyg != "homozygous_ref"
            rows.append({
                "gene": v.get("gene", gene_name),
                "rsid": v.get("rsid", ""),
                "variant_name": v.get("variant_name", ""),
                "genotype": v.get("genotype", ""),
                "zygosity": zyg,
                "clinical_significance": v.get("clinical_significance", ""),
                "drugs": v.get("drugs", ""),
                "evidence_level": v.get("evidence_level", ""),
                "is_variant": is_variant,
                "source": "Extended Panel",
            })
    return rows


def _collect_custom_pgx_results(pgx_dir: str) -> List[Dict[str, Any]]:
    """Read ``pgx_custom_result.json`` from the analysis pgx/ directory."""
    data = _read_pgx_custom_json(pgx_dir)
    if not data:
        return []
    return _custom_gene_rows_from_pgx_dict(data)


APOE_TAG_RSIDS = frozenset({"rs429358", "rs7412"})


def _zygosity_is_heterozygous(row: Dict[str, Any]) -> Optional[bool]:
    """Return True if het, False if homozygous (ref or alt), None if unknown."""
    z = (row.get("zygosity") or "").strip().lower()
    if z in ("heterozygous", "het"):
        return True
    if z in ("homozygous_ref", "homozygous_alt", "hom_ref", "hom_alt"):
        return False
    gt = (row.get("genotype") or "").strip().upper()
    if not gt:
        return None
    for sep in ("/", "|"):
        if sep in gt:
            parts = [p.strip() for p in gt.split(sep) if p.strip()]
            if len(parts) >= 2:
                return parts[0] != parts[1]
    return None


def _pipeline_claims_phased_apoe(raw: Optional[Dict[str, Any]]) -> Tuple[bool, str]:
    """True if upstream JSON explicitly reports phased APOE / haplotype (read-backed, trio, etc.)."""
    if not isinstance(raw, dict):
        return False, ""
    block = raw.get("apoe_phasing")
    if isinstance(block, dict):
        if block.get("phase_resolved") is True or block.get("phased") is True:
            hap = block.get("haplotype") or block.get("diplotype") or block.get("method")
            return True, str(hap or "pipeline")
    genes = raw.get("genes")
    if isinstance(genes, dict):
        for v in genes.get("APOE") or []:
            if isinstance(v, dict) and (v.get("phased") is True or v.get("phase_resolved") is True):
                return True, str(v.get("haplotype") or v.get("note") or "per-row")
    return False, ""


def apoe_phasing_assessment(
    custom_rows: List[Dict[str, Any]],
    raw_custom_json: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    APOE ε2/ε3/ε4 tag SNPs (rs429358, rs7412) — cis/trans cannot be resolved from **unphased**
    short-read exome when **both** loci are heterozygous. Set ``show_alert`` so UI/PDF can warn.

    If ``pgx_custom_result.json`` includes ``apoe_phasing`` with ``phased`` / ``phase_resolved``,
    or per-row ``phased`` on APOE variants, we treat phase as supplied by the pipeline.
    """
    apoe_rows = [
        r
        for r in custom_rows
        if isinstance(r, dict)
        and str(r.get("gene") or "").strip().upper() == "APOE"
        and str(r.get("rsid") or "").strip() in APOE_TAG_RSIDS
    ]
    by_rs = {str(r.get("rsid") or "").strip(): r for r in apoe_rows}

    phased, src = _pipeline_claims_phased_apoe(raw_custom_json)
    if phased:
        return {
            "status": "pipeline_resolved",
            "short_warning": "",
            "detail": (
                f"APOE phasing is marked as resolved by the pipeline ({src}). "
                "Confirm methodology in pgx_custom_result.json / lab SOP."
            ),
            "show_alert": False,
            "pipeline_phased": True,
        }

    if not apoe_rows:
        return {
            "status": "not_applicable",
            "short_warning": "",
            "detail": "",
            "show_alert": False,
        }

    r358 = by_rs.get("rs429358")
    r412 = by_rs.get("rs7412")
    if not r358 or not r412:
        return {
            "status": "incomplete",
            "short_warning": (
                "APOE: only one of the two tag SNPs (rs429358, rs7412) is present in extended PGx output."
            ),
            "detail": (
                "Full ε2/ε3/ε4 context usually requires both loci in pgx_custom_result.json. "
                "Verify the pipeline emitted both rows."
            ),
            "show_alert": True,
        }

    h358 = _zygosity_is_heterozygous(r358)
    h412 = _zygosity_is_heterozygous(r412)
    if h358 is True and h412 is True:
        return {
            "status": "ambiguous",
            "short_warning": (
                "APOE: rs429358 and rs7412 are both heterozygous — ε2/ε3/ε4 haplotype phase "
                "(cis vs trans) cannot be determined from unphased exome reads alone."
            ),
            "detail": (
                "Short-read whole-exome data does not resolve which alleles sit on the same chromosome "
                "when both tag SNPs are heterozygous. Do not report a definitive ε2/ε3/ε4 diplotype from "
                "this pattern alone without read-backed phasing, trio/family data, or orthogonal typing. "
                "Follow your laboratory’s policy for pharmacogenomic vs neurodegenerative risk reporting."
            ),
            "show_alert": True,
            "rs429358_heterozygous": True,
            "rs7412_heterozygous": True,
        }

    if h358 is None or h412 is None:
        return {
            "status": "unknown_zygosity",
            "short_warning": (
                "APOE: zygosity could not be determined for rs429358 and/or rs7412; phase assessment is limited."
            ),
            "detail": "Expected zygosity or genotype (e.g. 0/1, C/T) in pgx_custom_result.json.",
            "show_alert": True,
        }

    return {
        "status": "likely_unambiguous",
        "short_warning": "",
        "detail": (
            "At least one tag SNP is homozygous (or zygosity is not heterozygous at both loci), so "
            "ε2/ε3/ε4 diplotype is often inferable without cross-SNP phasing; still confirm per clinical standards."
        ),
        "show_alert": False,
        "rs429358_heterozygous": h358,
        "rs7412_heterozygous": h412,
    }


# ── APOE proactive PDF (ε2/ε3/ε4) — HTML fragments for customer report ─────────

_EPS = ("ε2", "ε3", "ε4")
_EPS_IDX = {e: i for i, e in enumerate(_EPS)}


def _sort_epsilon_diplotype(a: str, b: str) -> str:
    aa, bb = a.strip(), b.strip()
    if aa in _EPS_IDX and bb in _EPS_IDX:
        oa, ob = sorted([aa, bb], key=lambda x: _EPS_IDX[x])
        return f"{oa}/{ob}"
    return f"{aa}/{bb}"


def _normalize_pipeline_diplotype_key(s: str) -> Optional[str]:
    """Map pipeline free-text (e.g. 'ε3/ε4', '3/4') to canonical report_key."""
    if not s or not isinstance(s, str):
        return None
    t = s.strip()
    t = t.replace("ε", "").replace("e", "").replace("E", "")
    m = re.search(r"([234])\s*[/\s]\s*([234])", t)
    if not m:
        return None
    em = {"2": "ε2", "3": "ε3", "4": "ε4"}
    a, b = em.get(m.group(1)), em.get(m.group(2))
    if not a or not b:
        return None
    return _sort_epsilon_diplotype(a, b)


def _pipeline_explicit_diplotype(raw: Optional[Dict[str, Any]]) -> Optional[str]:
    if not isinstance(raw, dict):
        return None
    block = raw.get("apoe_phasing")
    if isinstance(block, dict):
        for k in ("diplotype", "haplotype", "epsilon_diplotype"):
            v = block.get(k)
            if v:
                nk = _normalize_pipeline_diplotype_key(str(v))
                if nk:
                    return nk
    genes = raw.get("genes")
    if isinstance(genes, dict):
        for v in genes.get("APOE") or []:
            if not isinstance(v, dict):
                continue
            for k in ("diplotype", "haplotype", "epsilon_diplotype"):
                val = v.get(k)
                if val:
                    nk = _normalize_pipeline_diplotype_key(str(val))
                    if nk:
                        return nk
    return None


def _parse_diploid_gt(gt: str) -> Tuple[str, str]:
    if not gt or not isinstance(gt, str):
        return "", ""
    s = gt.strip().upper().replace(" ", "")
    for sep in ("/", "|", ","):
        if sep in s:
            parts = s.split(sep, 1)
            if len(parts) >= 2:
                a, b = parts[0].strip(), parts[1].strip()
                ca = a[-1] if a else ""
                cb = b[-1] if b else ""
                if ca in "ACGT" and cb in "ACGT":
                    return ca, cb
    if len(s) >= 2 and s[0] in "ACGT" and s[1] in "ACGT":
        return s[0], s[1]
    return "", ""


def _haplo_epsilon(a358: str, a412: str) -> Optional[str]:
    """Map (rs429358 allele, rs7412 allele) to ε2/ε3/ε4 (GRCh38 tag-SNP convention)."""
    x = (a358 or "").upper()
    y = (a412 or "").upper()
    if x == "T" and y == "C":
        return "ε2"
    if x == "T" and y == "T":
        return "ε3"
    if x == "C" and y == "T":
        return "ε4"
    return None


def _zygosity_hint_heterozygous(row: Optional[Dict[str, Any]]) -> Optional[bool]:
    if not isinstance(row, dict):
        return None
    z = (row.get("zygosity") or "").strip().lower()
    if z in ("heterozygous", "het"):
        return True
    if z in ("homozygous_ref", "homozygous_alt", "hom_ref", "hom_alt"):
        return False
    g1, g2 = _parse_diploid_gt(row.get("genotype") or "")
    if g1 and g2:
        return g1 != g2
    return None


def infer_apoe_diplotype_for_report(
    custom_rows: List[Dict[str, Any]],
    raw_custom: Optional[Dict[str, Any]],
    apoe_phasing: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Canonical ε2/ε3/ε4 diplotype key for proactive PDF text, or ambiguous/unknown.

    ``report_key`` matches :data:`APOE_PROACTIVE_DIPLOTYPE_BODIES` (plus ``ambiguous_both_het``, ``unknown``).
    """
    pipe = _pipeline_explicit_diplotype(raw_custom)
    if pipe:
        return {"report_key": pipe, "source": "pipeline"}

    apoe_rows = [
        r
        for r in custom_rows
        if isinstance(r, dict)
        and str(r.get("gene") or "").strip().upper() == "APOE"
        and str(r.get("rsid") or "").strip() in APOE_TAG_RSIDS
    ]
    by_rs = {str(r.get("rsid") or "").strip(): r for r in apoe_rows}
    r358 = by_rs.get("rs429358")
    r412 = by_rs.get("rs7412")
    if not r358 or not r412:
        return {"report_key": "unknown", "source": "inferred", "reason": "incomplete_snps"}

    g358 = _parse_diploid_gt(str(r358.get("genotype") or ""))
    g412 = _parse_diploid_gt(str(r412.get("genotype") or ""))
    if (not g358[0] or not g412[0]) and isinstance(apoe_phasing, dict):
        if apoe_phasing.get("status") == "ambiguous":
            return {"report_key": "ambiguous_both_het", "source": "inferred"}
        if apoe_phasing.get("status") == "unknown_zygosity":
            return {"report_key": "unknown", "source": "inferred"}

    if not g358[0] or not g412[0]:
        # Fall back to zygosity-only inference
        h358 = _zygosity_hint_heterozygous(r358)
        h412 = _zygosity_hint_heterozygous(r412)
        if h358 is True and h412 is True:
            return {"report_key": "ambiguous_both_het", "source": "inferred"}
        return {"report_key": "unknown", "source": "inferred", "reason": "missing_genotype"}

    a358, b358 = g358
    a412, b412 = g412
    c358 = Counter([a358, b358])
    c412 = Counter([a412, b412])
    u358 = set(c358.keys())
    u412 = set(c412.keys())

    if len(u358) == 1 and len(u412) == 1:
        h = _haplo_epsilon(a358, a412)
        if h:
            return {"report_key": f"{h}/{h}", "source": "inferred"}
        return {"report_key": "unknown", "source": "inferred"}

    # Both loci heterozygous (unphased): chromosomes may pair as ε2+ε4 or ε3+ε3 — not resolved from exome alone.
    if len(u358) == 2 and len(u412) == 2:
        return {"report_key": "ambiguous_both_het", "source": "inferred"}

    if len(u358) == 2 and len(u412) == 1:
        al412 = next(iter(u412))
        h1 = _haplo_epsilon(a358, al412)
        h2 = _haplo_epsilon(b358, al412)
        if h1 and h2:
            return {"report_key": _sort_epsilon_diplotype(h1, h2), "source": "inferred"}
        return {"report_key": "unknown", "source": "inferred"}

    if len(u358) == 1 and len(u412) == 2:
        al358 = next(iter(u358))
        h1 = _haplo_epsilon(al358, a412)
        h2 = _haplo_epsilon(al358, b412)
        if h1 and h2:
            return {"report_key": _sort_epsilon_diplotype(h1, h2), "source": "inferred"}
        return {"report_key": "unknown", "source": "inferred"}

    return {"report_key": "unknown", "source": "inferred"}


APOE_PROACTIVE_DISCLAIMER_HTML = (
    '<p style="margin-top:12px;padding-top:10px;border-top:1px solid #cbd5e1;font-size:7.5pt;line-height:1.4;color:#475569">'
    "<strong>Disclaimer:</strong> Risk estimates are approximate and based on population studies. "
    "Actual risk varies depending on age, sex, ancestry, environmental factors, and family history. "
    "APOE genotype is not diagnostic and should not be used alone to predict disease."
    "</p>"
)

# Single-diplotype bodies only (paired title + risk + clinical); disclaimer appended separately.
APOE_PROACTIVE_DIPLOTYPE_BODIES: Dict[str, str] = {
    "ε2/ε2": (
        "<p style=\"margin:0 0 6px;font-size:11pt;font-weight:700\">ε2 / ε2</p>"
        '<p style="margin:0 0 4px"><strong>Alzheimer disease risk:</strong><br />~0.5× (reduced risk)</p>'
        '<p style="margin:0"><strong>Clinical summary:</strong><br />'
        "Associated with reduced risk of Alzheimer disease. May be associated with type III "
        "hyperlipoproteinemia in some individuals.</p>"
    ),
    "ε2/ε3": (
        "<p style=\"margin:0 0 6px;font-size:11pt;font-weight:700\">ε2 / ε3</p>"
        '<p style="margin:0 0 4px"><strong>Alzheimer disease risk:</strong><br />~0.6–0.8× (slightly reduced)</p>'
        '<p style="margin:0"><strong>Clinical summary:</strong><br />'
        "Generally considered protective or neutral. Mild effects on lipid metabolism may be observed.</p>"
    ),
    "ε3/ε3": (
        "<p style=\"margin:0 0 6px;font-size:11pt;font-weight:700\">ε3 / ε3</p>"
        '<p style="margin:0 0 4px"><strong>Alzheimer disease risk:</strong><br />~1× (baseline)</p>'
        '<p style="margin:0"><strong>Clinical summary:</strong><br />'
        "Represents the reference genotype with average population risk.</p>"
    ),
    "ε2/ε4": (
        "<p style=\"margin:0 0 6px;font-size:11pt;font-weight:700\">ε2 / ε4</p>"
        '<p style="margin:0 0 4px"><strong>Alzheimer disease risk:</strong><br />'
        "<strong>Intermediate</strong> (approximate relative risk often cited ~2–3× vs ε3/ε3 baseline; population estimates vary)</p>"
        '<p style="margin:0"><strong>Clinical summary:</strong><br />'
        "Mixed genotype with variable risk. ε4 increases risk, ε2 may partially offset.</p>"
    ),
    "ambiguous_both_het": (
        '<p style="margin:0 0 8px;font-size:10.5pt;font-weight:700">APOE haplotype — phase ambiguity</p>'
        '<p style="margin:0 0 8px;font-size:8.5pt;line-height:1.45">'
        "<strong>APOE haplotype could not be definitively determined due to phase ambiguity.</strong></p>"
        '<p style="margin:0 0 6px;font-size:8.5pt;line-height:1.45">'
        "The detected variants at the tag SNPs (rs429358, rs7412) are consistent with <strong>either</strong>:</p>"
        '<ul style="margin:0 0 10px;padding-left:18px;font-size:8.5pt;line-height:1.45">'
        "<li>ε2/ε4 genotype, <em>or</em></li>"
        "<li>ε3/ε3 genotype.</li>"
        "</ul>"
        '<p style="margin:0;font-size:8.5pt;line-height:1.45;color:#334155">'
        "Additional testing (e.g., targeted genotyping or long-read sequencing) may be considered to "
        "resolve haplotype phase if clinically indicated.</p>"
    ),
    "ε3/ε4": (
        "<p style=\"margin:0 0 6px;font-size:11pt;font-weight:700\">ε3 / ε4</p>"
        '<p style="margin:0 0 4px"><strong>Alzheimer disease risk:</strong><br />'
        "<strong>Moderate increase</strong> (approximate relative risk often cited ~2–4× vs ε3/ε3; population estimates vary)</p>"
        '<p style="margin:0"><strong>Clinical summary:</strong><br />'
        "Associated with moderately increased risk of Alzheimer disease. May also be associated with "
        "higher LDL cholesterol.</p>"
    ),
    "ε4/ε4": (
        "<p style=\"margin:0 0 6px;font-size:11pt;font-weight:700\">ε4 / ε4</p>"
        '<p style="margin:0 0 4px"><strong>Alzheimer disease risk:</strong><br />'
        "<strong>High risk</strong> (approximate relative risk often cited ~8–12× vs ε3/ε3; population estimates vary)</p>"
        '<p style="margin:0"><strong>Clinical summary:</strong><br />'
        "Associated with significantly increased risk and earlier onset. Also linked to increased "
        "cardiovascular risk.</p>"
    ),
    "unknown": (
        '<p style="margin:0 0 6px;font-size:10pt"><strong>APOE ε2/ε3/ε4 summary</strong></p>'
        '<p style="margin:0;font-size:8.5pt">Diplotype could not be determined from tag-SNP genotypes '
        "in this report. Refer to raw PGx outputs and laboratory SOP.</p>"
    ),
}


def _apoe_proactive_pdf_alert_html(report_key: str) -> str:
    """Prominent alert strip for customer PDF (risk tier or unresolved phase)."""
    rk = (report_key or "").strip()
    if rk in ("ε2/ε4", "ε3/ε4", "ε4/ε4"):
        labels = {
            "ε2/ε4": ("Intermediate", "#fffbeb", "#d97706", "#92400e"),
            "ε3/ε4": ("Moderate increase", "#fffbeb", "#d97706", "#92400e"),
            "ε4/ε4": ("High risk", "#fef2f2", "#dc2626", "#991b1b"),
        }
        label, bg, border, fg = labels[rk]
        disp = rk.replace("/", " / ")
        return (
            f'<div style="margin:0 0 12px;padding:10px 12px;border-radius:8px;border:2px solid {border};'
            f'background:{bg};font-size:9pt;line-height:1.45;color:{fg}">'
            f"<strong>Clinical alert — APOE</strong><br />"
            f"Diplotype <strong>{html.escape(disp)}</strong>: Alzheimer disease risk — "
            f"<strong>{html.escape(label)}</strong>. Interpret in full clinical context.</div>"
        )
    if rk in ("ambiguous_both_het", "unknown"):
        return (
            '<div style="margin:0 0 12px;padding:10px 12px;border-radius:8px;border:2px solid #dc2626;'
            'background:#fef2f2;font-size:9pt;line-height:1.45;color:#991b1b">'
            "<strong>Phasing alert — APOE</strong><br />"
            "ε2/ε3/ε4 haplotype phase could not be determined from the available tag-SNP data. "
            "Do not report a definitive diplotype for preventive-health counseling without resolving phase "
            "per laboratory policy (e.g. orthogonal typing or read-backed phasing).</div>"
        )
    return ""


def build_apoe_proactive_pdf_html(report_key: str) -> str:
    """WeasyPrint-safe HTML fragment (optional alert + one diplotype block + disclaimer)."""
    rk = (report_key or "").strip() or "unknown"
    alert = _apoe_proactive_pdf_alert_html(rk)
    body = APOE_PROACTIVE_DIPLOTYPE_BODIES.get(rk) or APOE_PROACTIVE_DIPLOTYPE_BODIES["unknown"]
    return (
        '<div class="apoe-proactive-pdf" style="font-size:8.5pt;line-height:1.45;color:#0f172a">'
        f"{alert}{body}{APOE_PROACTIVE_DISCLAIMER_HTML}</div>"
    )


def load_pgx_custom_variants_reference() -> List[Dict[str, str]]:
    """Load the built-in pgx_custom_variants.tsv reference file."""
    ref_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "data", "db", "pgx_custom_variants.tsv",
    )
    ref_path = os.path.normpath(ref_path)
    if not os.path.isfile(ref_path):
        return []
    rows: List[Dict[str, str]] = []
    try:
        with open(ref_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 8:
                    continue
                rows.append({
                    "gene": cols[0],
                    "rsid": cols[1],
                    "variant_name": cols[6],
                    "clinical_significance": cols[7],
                    "drugs": cols[8] if len(cols) > 8 else "",
                    "evidence_level": cols[9] if len(cols) > 9 else "",
                    "clinpgx_url": cols[10] if len(cols) > 10 else "",
                })
    except Exception as e:
        logger.warning("[pgx] could not load pgx_custom_variants.tsv: %s", e)
    return rows


def filter_pgx_gene_results_by_panel(
    gene_results: List[Dict[str, Any]],
    panel_genes: AbstractSet[str],
) -> List[Dict[str, Any]]:
    """
    Keep only PGx rows whose ``gene`` is in the order interpretation set (WES panel + extras).

    Same gene universe as variant post-filter: ``interpretation_gene_set_for_job`` / ``panel_interpretation_genes``.
    When ``panel_genes`` is empty, returns ``gene_results`` unchanged (caller should not pass an empty set to mean “filter”).
    """
    if not panel_genes:
        return list(gene_results)
    up = {str(g).strip().upper() for g in panel_genes if g is not None and str(g).strip()}
    if not up:
        return list(gene_results)
    out: List[Dict[str, Any]] = []
    for row in gene_results:
        if not isinstance(row, dict):
            continue
        g = row.get("gene")
        if g is None or str(g).strip().upper() not in up:
            continue
        out.append(row)
    return out


def merge_pgx_gene_reviews(
    fresh_rows: List[Dict[str, Any]],
    previous_pgx: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """After reprocess, keep ``reviewer_confirmed`` / ``reviewer_comment`` per ``gene``."""
    prev_by_gene: Dict[str, Dict[str, Any]] = {}
    if isinstance(previous_pgx, dict):
        for row in previous_pgx.get("gene_results") or []:
            if isinstance(row, dict) and row.get("gene"):
                prev_by_gene[str(row["gene"])] = row
    out: List[Dict[str, Any]] = []
    for row in fresh_rows:
        if not isinstance(row, dict):
            continue
        g = row.get("gene")
        merged = dict(row)
        pr = prev_by_gene.get(str(g)) if g is not None else None
        if pr:
            merged["reviewer_confirmed"] = bool(pr.get("reviewer_confirmed", False))
            merged["reviewer_comment"] = (pr.get("reviewer_comment") or "")[:4000]
        else:
            merged.setdefault("reviewer_confirmed", False)
            merged.setdefault("reviewer_comment", "")
        out.append(merged)
    return out


def merge_pgx_custom_gene_reviews(
    fresh_rows: List[Dict[str, Any]],
    previous_pgx: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """After reprocess, keep ``reviewer_confirmed`` / ``reviewer_comment`` per ``gene`` + ``rsid`` (extended panel / APOE)."""
    prev_by_key: Dict[Tuple[str, str], Dict[str, Any]] = {}
    if isinstance(previous_pgx, dict):
        for row in previous_pgx.get("custom_gene_results") or []:
            if isinstance(row, dict):
                g = str(row.get("gene") or "").strip()
                rs = str(row.get("rsid") or "").strip()
                if g and rs:
                    prev_by_key[(g, rs)] = row
    out: List[Dict[str, Any]] = []
    for row in fresh_rows:
        if not isinstance(row, dict):
            continue
        g = str(row.get("gene") or "").strip()
        rs = str(row.get("rsid") or "").strip()
        merged = dict(row)
        pr = prev_by_key.get((g, rs)) if g and rs else None
        if pr:
            merged["reviewer_confirmed"] = bool(pr.get("reviewer_confirmed", False))
            merged["reviewer_comment"] = (pr.get("reviewer_comment") or "")[:4000]
        else:
            merged.setdefault("reviewer_confirmed", False)
            merged.setdefault("reviewer_comment", "")
        out.append(merged)
    return out


def merge_pgx_portal_review(
    fresh: Dict[str, Any],
    previous: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    When ``generate_result_json`` rewrites ``pgx``, keep ``portal_review`` from the prior
    ``result.json`` (reviewer notes / reviewed flag from the portal).
    """
    out = dict(fresh) if isinstance(fresh, dict) else {}
    if not isinstance(previous, dict):
        return out
    pr = previous.get("portal_review")
    if isinstance(pr, dict) and pr:
        out["portal_review"] = dict(pr)
    ap = previous.get("apoe_phasing")
    if isinstance(ap, dict) and not out.get("apoe_phasing"):
        out["apoe_phasing"] = dict(ap)
    return out


def collect_pgx_from_analysis_dir(root: str, sample_name: str) -> Dict[str, Any]:
    """
    Scan ``root/pgx/`` for PharmCAT meta + summary text.

    Returns a dict suitable for ``result.json["pgx"]``.
    """
    root = (root or "").strip()
    if not root or not os.path.isdir(root):
        return {
            "status": "not_found",
            "message": "analysis/output root is missing or not a directory",
            "searched": root,
        }

    pgx_dir = os.path.join(root, "pgx")
    if not os.path.isdir(pgx_dir):
        return {
            "status": "not_found",
            "message": f"No pgx directory under {root}",
            "searched_root": os.path.abspath(root),
        }

    meta: Dict[str, Any] = {}
    meta_path = os.path.join(pgx_dir, "pgx_meta.json")
    if os.path.isfile(meta_path):
        try:
            with open(meta_path, "r", encoding="utf-8") as f:
                meta = json.load(f)
        except Exception as e:
            logger.warning("[pgx] could not read %s: %s", meta_path, e)

    summary_text = ""
    summary_path = os.path.join(pgx_dir, "pgx_summary.txt")
    if os.path.isfile(summary_path):
        try:
            with open(summary_path, "r", encoding="utf-8") as f:
                summary_text = f.read()
        except Exception as e:
            logger.warning("[pgx] could not read %s: %s", summary_path, e)

    reporter_html = _first_basename(pgx_dir, "*_pgx.report.html")
    result_json_name: Optional[str] = None
    rp = os.path.join(pgx_dir, "pgx_result.json")
    if os.path.isfile(rp):
        result_json_name = "pgx_result.json"

    exit_ok = str((meta or {}).get("exit_status") or "").lower() in ("success", "ok", "0")
    if not summary_text.strip():
        if exit_ok and meta:
            summary_text = (
                f"PGx completed ({(meta.get('tool_version') or 'PharmCAT')}). "
                "Summary file was empty; see full PharmCAT HTML in pipeline pgx/ output."
            )
        else:
            return {
                "status": "error",
                "message": "pgx_summary.txt missing or empty",
                "meta": meta,
                "pgx_dir": os.path.abspath(pgx_dir),
            }

    out: Dict[str, Any] = {
        "status": "ok",
        "summary_text": summary_text,
        "meta": meta,
        "artifacts": {
            "pgx_summary_txt": "pgx/pgx_summary.txt",
            "pgx_meta_json": "pgx/pgx_meta.json",
            "reporter_html_basename": reporter_html,
            "pgx_result_json": result_json_name,
        },
        "pgx_dir": os.path.abspath(pgx_dir),
        "gene_results": [],
    }

    rp_full = os.path.join(pgx_dir, "pgx_result.json")
    if os.path.isfile(rp_full):
        try:
            with open(rp_full, "r", encoding="utf-8") as f:
                bundle = json.load(f)
            phen = bundle.get("phenotype") if isinstance(bundle, dict) else None
            if isinstance(phen, dict):
                out["gene_results"] = extract_gene_results_from_phenotype(phen)
                gr = phen.get("geneReports")
                if isinstance(gr, dict):
                    if "CPIC" in gr or "DPWG" in gr:
                        pharmcat_all = set()
                        for src in ("CPIC", "DPWG"):
                            rr = gr.get(src)
                            if isinstance(rr, dict):
                                pharmcat_all |= set(rr.keys())
                    else:
                        pharmcat_all = {k for k, v in gr.items() if isinstance(v, dict)}
                    out["all_pharmcat_genes"] = sorted(pharmcat_all)
        except Exception as e:
            logger.warning("[pgx] could not parse gene rows from %s: %s", rp_full, e)

    actionable_genes: set = {
        r["gene"] for r in out.get("gene_results", [])
        if isinstance(r, dict) and r.get("category") == "actionable" and r.get("gene")
    }
    if actionable_genes:
        drug_recs, drug_src = _collect_drug_recommendations_from_pgx_dir(
            pgx_dir, actionable_genes
        )
        out["drug_recommendations"] = drug_recs
        if drug_src:
            out["artifacts"]["drug_recommendations_source"] = drug_src
    else:
        out["drug_recommendations"] = []

    raw_custom = _read_pgx_custom_json(pgx_dir)
    custom = _custom_gene_rows_from_pgx_dict(raw_custom) if raw_custom else []
    if custom:
        out["custom_gene_results"] = custom
    try:
        out["apoe_phasing"] = apoe_phasing_assessment(custom, raw_custom)
        out["apoe_diplotype_for_report"] = infer_apoe_diplotype_for_report(
            custom, raw_custom, out.get("apoe_phasing")
        )
    except Exception as e:
        logger.warning("[pgx] apoe phasing/diplotype step failed (non-fatal): %s", e)
        out.setdefault("apoe_phasing", {"status": "error", "show_alert": False})
        out["apoe_diplotype_for_report"] = {
            "report_key": "unknown",
            "source": "error",
            "message": str(e),
        }

    return out


def _first_basename(dir_path: str, pattern: str) -> Optional[str]:
    paths = sorted(glob.glob(os.path.join(dir_path, pattern)))
    if not paths:
        return None
    return os.path.basename(paths[0])


def pgx_for_pdf(pgx: Dict[str, Any]) -> Dict[str, Any]:
    """
    Subset + HTML fragment for WeasyPrint (written into report.json used by Jinja).
    """
    if not isinstance(pgx, dict):
        return {}
    st = pgx.get("status")
    if st == "not_found":
        return {}

    meta = pgx.get("meta") if isinstance(pgx.get("meta"), dict) else {}
    tool_v = (meta.get("tool_version") or meta.get("tool") or "PharmCAT").strip()

    if st == "error":
        msg = html.escape(str(pgx.get("message") or "PGx error"))
        err_out: Dict[str, Any] = {
            **pgx,
            "summary_for_pdf_html": (
                f'<p class="detail-text" style="color:#b91c1c;font-size:9pt;">{msg}</p>'
            ),
            "tool_version_line": tool_v,
        }
        err_out.pop("portal_review", None)
        err_out.pop("gene_results", None)
        return err_out

    body = pgx.get("summary_text") or ""
    esc = html.escape(body)
    foot = (
        f'<p class="muted" style="font-size:7.5pt;margin-top:10px;line-height:1.4;">'
        f"Source: PharmCAT summary text. Full interactive report: "
        f"<code>pgx/</code> output (<code>*_pgx.report.html</code>). {html.escape(tool_v)}"
        f"</p>"
    )
    apoe_html = ""
    apoe_ph = pgx.get("apoe_phasing")
    if isinstance(apoe_ph, dict) and apoe_ph.get("show_alert"):
        sw = html.escape(str(apoe_ph.get("short_warning") or "").strip())
        det = html.escape(str(apoe_ph.get("detail") or "").strip())
        if sw or det:
            br_sw = f"<br />{sw}" if sw else ""
            br_det = (
                f'<br /><span style="font-size:8pt;opacity:.95">{det}</span>' if det else ""
            )
            apoe_html = (
                f'<div class="pgx-apoe-phase" style="margin:0 0 12px;padding:10px 12px;border-radius:8px;'
                f"border:1px solid #f59e0b;background:#fffbeb;font-size:8.5pt;line-height:1.45;color:#92400e\">"
                f"<strong>APOE phasing</strong>"
                f"{br_sw}{br_det}</div>"
            )
    block = (
        f"{apoe_html}"
        f'<pre class="pgx-summary" style="white-space:pre-wrap;font-family:ui-monospace,monospace;'
        f"font-size:8.5pt;line-height:1.35;border:1px solid #cbd5e1;border-radius:8px;"
        f'padding:12px;background:#f8fafc;">{esc}</pre>{foot}'
    )
    out = dict(pgx)
    out["summary_for_pdf_html"] = block
    out["tool_version_line"] = tool_v

    pr_pdf = pgx.get("portal_review") if isinstance(pgx.get("portal_review"), dict) else {}

    def _apoe_tag_rows_present() -> bool:
        cgr = pgx.get("custom_gene_results")
        if not isinstance(cgr, list):
            return False
        tag = {"rs429358", "rs7412"}
        for r in cgr:
            if not isinstance(r, dict):
                continue
            if str(r.get("gene") or "").strip().upper() != "APOE":
                continue
            if str(r.get("rsid") or "").strip() in tag:
                return True
        return False

    ex_apoe = pr_pdf.get("include_apoe_proactive_pdf")
    if ex_apoe is False or str(ex_apoe).strip().lower() in ("false", "0", "no"):
        include_apoe_pdf = False
    elif ex_apoe is True or str(ex_apoe).strip().lower() in ("true", "1", "yes"):
        include_apoe_pdf = True
    else:
        # Reviewer never saved PGx review: still emit APOE PDF when extended panel has tag SNPs
        include_apoe_pdf = _apoe_tag_rows_present()

    adp = pgx.get("apoe_diplotype_for_report")
    if include_apoe_pdf and not isinstance(adp, dict):
        cgr = pgx.get("custom_gene_results")
        if isinstance(cgr, list):
            adp = infer_apoe_diplotype_for_report(
                cgr, None, pgx.get("apoe_phasing") if isinstance(pgx.get("apoe_phasing"), dict) else None
            )
    if include_apoe_pdf and isinstance(adp, dict):
        rk = str(adp.get("report_key") or "unknown").strip()
        out["apoe_proactive_summary_html"] = build_apoe_proactive_pdf_html(rk)
    else:
        out["apoe_proactive_summary_html"] = ""

    out.pop("portal_review", None)

    # Build the full gene list from ALL sources BEFORE reviewer_confirmed filtering
    all_gene_names: set = set()
    all_genes = out.get("gene_results") if isinstance(out.get("gene_results"), list) else []
    for r in all_genes:
        if isinstance(r, dict) and r.get("gene"):
            all_gene_names.add(str(r["gene"]).strip())
    pharmcat_all = pgx.get("all_pharmcat_genes")
    if isinstance(pharmcat_all, list):
        for g in pharmcat_all:
            if g:
                all_gene_names.add(str(g).strip())
    custom = pgx.get("custom_gene_results")
    if isinstance(custom, list):
        for r in custom:
            if isinstance(r, dict) and r.get("gene"):
                all_gene_names.add(str(r["gene"]).strip())
    _EXCLUDE_FROM_GENE_LIST = {"MT-RNR1", "CFTR", "HLA-A", "HLA-B", "CES1", "IFNL3"}
    out["genes_evaluated"] = sorted(all_gene_names - _EXCLUDE_FROM_GENE_LIST)

    confirmed = [r for r in all_genes if isinstance(r, dict) and r.get("reviewer_confirmed")]
    out["gene_results"] = confirmed if confirmed else all_genes

    drug_recs = pgx.get("drug_recommendations")
    if isinstance(drug_recs, list) and drug_recs:
        confirmed_gene_names = {r["gene"] for r in out["gene_results"] if r.get("gene")}
        out["drug_recommendations"] = [
            r for r in drug_recs
            if not confirmed_gene_names or r.get("gene") in confirmed_gene_names
        ]
    else:
        out["drug_recommendations"] = generate_cpic_drug_recommendations(out["gene_results"])

    if isinstance(custom, list) and custom:
        pharmcat_genes = {r["gene"] for r in out["gene_results"] if r.get("gene")}
        non_pharmcat = [r for r in custom if r.get("gene") not in pharmcat_genes]
        confirmed_custom = [r for r in non_pharmcat if r.get("reviewer_confirmed")]
        # Match PharmCAT row logic: if reviewer checked at least one ✓ Include, PDF shows only those;
        # if none checked, show all extended-panel rows so the PGx report is not empty.
        out["custom_gene_results"] = confirmed_custom if confirmed_custom else non_pharmcat
    else:
        out["custom_gene_results"] = []

    return out


def sanitize_pgx_payload_for_pdf_render(report_data: Dict[str, Any]) -> None:
    """Ensure ``report_data['pgx']`` has ``summary_for_pdf_html`` when raw summary exists."""
    pgx = report_data.get("pgx")
    if not isinstance(pgx, dict):
        return
    if (pgx.get("summary_for_pdf_html") or "").strip():
        return
    merged = pgx_for_pdf(pgx)
    if merged:
        combined = {**pgx, **merged}
        combined.pop("portal_review", None)
        report_data["pgx"] = combined
