"""VEP CSQ transcript pick when SYMBOL empty on severity-first row (import vep_parser only)."""

import importlib.util
import os

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_vep_parser():
    path = os.path.join(_ROOT, "app", "services", "carrier_screening", "vep_parser.py")
    spec = importlib.util.spec_from_file_location("vep_parser_gene_pick", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


_vp = _load_vep_parser()


def test_parse_csq_prefers_row_with_symbol_when_top_is_empty():
    # Minimal CSQ field list; two transcripts same ALT "A", first MODIFIER + empty SYMBOL
    fields = [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",
        "BIOTYPE",
        "EXON",
        "INTRON",
        "HGVSc",
        "HGVSp",
    ]

    def row(sym, impact, cons, hgvsc=""):
        parts = [
            "A",
            cons,
            impact,
            sym,
            "ENSG00000000001",
            "Transcript",
            "ENST00000000001",
            "protein_coding",
            "",
            "",
            hgvsc,
            "",
        ]
        return "|".join(parts)

    csq = row("", "MODIFIER", "downstream_gene_variant") + "," + row(
        "TP53", "MODERATE", "missense_variant", "ENST00000000001:c.101A>G"
    )
    out = _vp.parse_csq_record(csq, fields, alt_allele="A", prefer_mane=False)
    assert out is not None
    assert out.get("gene") == "TP53"
    assert "c.101A>G" in (out.get("hgvsc") or "")


def test_normalize_falls_back_to_ensembl_gene_id():
    rec = {
        "SYMBOL": "",
        "Gene": "ENSG00000141510",
        "HGVSc": "",
        "HGVSp": "",
        "Consequence": "intron_variant",
        "IMPACT": "MODIFIER",
        "Feature": "",
        "BIOTYPE": "",
        "EXON": "",
        "INTRON": "",
        "gnomADe_AF": "",
        "gnomADg_AF": "",
        "Existing_variation": "",
        "SIFT": "",
        "PolyPhen": "",
        "CLIN_SIG": "",
        "MANE_SELECT": "",
        "CANONICAL": "",
    }
    out = _vp._normalize_csq_record(rec)
    assert out.get("gene") == "ENSG00000141510"
