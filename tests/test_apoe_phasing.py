"""APOE tag-SNP phasing hints for extended PGx (rs429358 / rs7412)."""

from app.services.carrier_screening.pgx_report import apoe_phasing_assessment


def test_apoe_both_heterozygous_is_ambiguous():
    rows = [
        {"gene": "APOE", "rsid": "rs429358", "zygosity": "heterozygous", "genotype": "T/C"},
        {"gene": "APOE", "rsid": "rs7412", "zygosity": "heterozygous", "genotype": "C/T"},
    ]
    r = apoe_phasing_assessment(rows, None)
    assert r["status"] == "ambiguous"
    assert r["show_alert"] is True


def test_pipeline_marks_phased():
    rows = [
        {"gene": "APOE", "rsid": "rs429358", "zygosity": "heterozygous"},
        {"gene": "APOE", "rsid": "rs7412", "zygosity": "heterozygous"},
    ]
    raw = {"apoe_phasing": {"phased": True, "haplotype": "ε3/ε4"}}
    r = apoe_phasing_assessment(rows, raw)
    assert r["status"] == "pipeline_resolved"
    assert r["show_alert"] is False


def test_one_homozygous_likely_unambiguous():
    rows = [
        {"gene": "APOE", "rsid": "rs429358", "zygosity": "homozygous_ref", "genotype": "T/T"},
        {"gene": "APOE", "rsid": "rs7412", "zygosity": "heterozygous", "genotype": "C/T"},
    ]
    r = apoe_phasing_assessment(rows, None)
    assert r["status"] == "likely_unambiguous"
    assert r["show_alert"] is False


def test_incomplete_loci():
    rows = [{"gene": "APOE", "rsid": "rs429358", "zygosity": "heterozygous"}]
    r = apoe_phasing_assessment(rows, None)
    assert r["status"] == "incomplete"
    assert r["show_alert"] is True
