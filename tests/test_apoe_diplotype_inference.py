"""APOE ε2/ε3/ε4 diplotype inference for proactive PDF text."""

from app.services.carrier_screening.pgx_report import infer_apoe_diplotype_for_report


def _apoe(r358_gt: str, r412_gt: str):
    return [
        {"gene": "APOE", "rsid": "rs429358", "genotype": r358_gt, "zygosity": "homozygous"},
        {"gene": "APOE", "rsid": "rs7412", "genotype": r412_gt, "zygosity": "homozygous"},
    ]


def test_infer_e3e3():
    r = infer_apoe_diplotype_for_report(_apoe("T/T", "T/T"), None, None)
    assert r["report_key"] == "ε3/ε3"
    assert r["source"] == "inferred"


def test_infer_e4e4():
    r = infer_apoe_diplotype_for_report(_apoe("C/C", "T/T"), None, None)
    assert r["report_key"] == "ε4/ε4"


def test_infer_e2e2():
    r = infer_apoe_diplotype_for_report(_apoe("T/T", "C/C"), None, None)
    assert r["report_key"] == "ε2/ε2"


def test_infer_e3e4_heterozygous_358_only():
    r = infer_apoe_diplotype_for_report(_apoe("T/C", "T/T"), None, None)
    assert r["report_key"] == "ε3/ε4"


def test_infer_e2e3_heterozygous_412_only():
    r = infer_apoe_diplotype_for_report(_apoe("T/T", "C/T"), None, None)
    assert r["report_key"] == "ε2/ε3"


def test_infer_ambiguous_both_heterozygous():
    r = infer_apoe_diplotype_for_report(_apoe("T/C", "C/T"), None, None)
    assert r["report_key"] == "ambiguous_both_het"


def test_pipeline_diplotype_overrides():
    raw = {"apoe_phasing": {"diplotype": "ε3/ε4"}}
    r = infer_apoe_diplotype_for_report([], raw, None)
    assert r["report_key"] == "ε3/ε4"
    assert r["source"] == "pipeline"
