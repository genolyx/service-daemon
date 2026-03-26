"""
Reference / annotation 파일 존재 여부 요약 (헬스·운영 점검용).
"""

from __future__ import annotations

import glob
import os
from typing import Any, Dict, List, Optional, Tuple

from .config import settings


def _stat_path(path: Optional[str], *, is_dir: bool = False) -> Dict[str, Any]:
    if not path or not str(path).strip():
        return {"path": None, "ok": False, "reason": "not configured"}
    p = os.path.abspath(str(path).strip())
    if is_dir:
        ok = os.path.isdir(p)
    else:
        ok = os.path.isfile(p)
    return {
        "path": p,
        "ok": ok,
        "reason": None if ok else ("not a directory" if is_dir else "not a file"),
    }


def annotation_resource_report() -> Dict[str, Any]:
    """
    compose 가 마운트하는 /data/reference 등과 .env 경로가 일치하는지 점검.
    """
    checks: List[Tuple[str, Optional[str], bool]] = [
        ("ref_fasta", settings.ref_fasta, False),
        ("ref_fai", settings.ref_fai, False),
        ("ref_dict", settings.ref_dict, False),
        ("ref_bwa_indices", settings.ref_bwa_indices, True),
        ("clinvar_vcf", settings.clinvar_vcf, False),
        ("gnomad_vcf", settings.gnomad_vcf, False),
        ("gnomad_dir", settings.gnomad_dir, True),
        ("dbsnp_vcf", settings.dbsnp_vcf, False),
        ("clingen_tsv", settings.clingen_tsv, False),
        ("mane_gff", settings.mane_gff, False),
        ("gene_bed", settings.gene_bed, False),
        ("hpo_gene_file", settings.hpo_gene_file, False),
        ("curated_variants_db", settings.curated_variants_db, False),
        ("hgmd_vcf", settings.hgmd_vcf, False),
        ("snpeff_jar", settings.snpeff_jar, False),
    ]

    resources: Dict[str, Any] = {}
    for key, path, is_dir in checks:
        resources[key] = _stat_path(path, is_dir=is_dir)

    data_dir = getattr(settings, "snpeff_data_dir", None) or ""
    if str(data_dir).strip():
        sd = _stat_path(data_dir, is_dir=True)
        resources["snpeff_data_dir"] = sd
        if sd["ok"]:
            gdir = os.path.join(sd["path"], "genomes")
            resources["snpeff_genomes_dir"] = _stat_path(gdir, is_dir=True)
    else:
        resources["snpeff_data_dir"] = {
            "path": None,
            "ok": False,
            "reason": "not configured (set SNPEFF_DATA_DIR for downloadable genomes)",
        }

    # gnomAD 디렉터리에 exome bgz 가 있는지 요약
    gdir = settings.gnomad_dir
    gnomad_files = []
    if gdir and os.path.isdir(gdir):
        for pat in (settings.gnomad_exomes_glob, settings.gnomad_genomes_glob):
            gnomad_files.extend(
                glob.glob(os.path.join(gdir, pat))[:5]
            )
    resources["gnomad_glob_sample_count"] = len(gnomad_files)

    dg = settings.disease_gene_json
    resources["disease_gene_json"] = _stat_path(dg, is_dir=False)

    if settings.carrier_default_backbone_bed:
        resources["carrier_default_backbone_bed"] = _stat_path(
            settings.carrier_default_backbone_bed, is_dir=False
        )
    if settings.carrier_default_disease_bed:
        resources["carrier_default_disease_bed"] = _stat_path(
            settings.carrier_default_disease_bed, is_dir=False
        )

    configured_but_missing = [
        k
        for k, v in resources.items()
        if k != "gnomad_glob_sample_count" and isinstance(v, dict) and v.get("path") and not v.get("ok")
    ]
    return {
        "resources": resources,
        "configured_but_missing_paths": configured_but_missing,
        "hint": (
            "Host: populate BASE_DIR/reference (see .env.docker comments). "
            "snpEff: install jar in image; run once: java -jar snpEff.jar download -v <DB> -dataDir $SNPEFF_DATA_DIR"
        ),
    }
