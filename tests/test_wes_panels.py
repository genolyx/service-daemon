"""WES/exome panel catalog (order-time BED selection)."""

from app.services.wes_panels import (
    apply_wes_panel_to_job_params,
    get_panel_by_id,
    interpretation_gene_set_for_job,
    list_panels,
    panels_for_api_response,
    should_apply_interpretation_post_filter,
    split_genes_from_text,
)
from app.services.carrier_screening.plugin import _filter_variants_by_interpretation_genes
from app.models import Job


def test_catalog_loads(wes_test_catalog):
    panels = list_panels()
    assert isinstance(panels, list)
    assert any(p.get("id") == "test_carrier_panel" for p in panels)
    # custom file may be empty; merged list still includes bundled
    assert any((p.get("_source") or "") in ("bundled", "custom") for p in panels)


def test_get_panel(wes_test_catalog):
    p = get_panel_by_id("test_carrier_panel")
    assert p is not None
    assert "disease_bed" in p


def test_apply_sets_disease_bed(wes_test_catalog):
    job = Job(
        order_id="x",
        service_code="carrier_screening",
        sample_name="x",
        work_dir="01",
        params={"wes_panel_id": "test_carrier_panel"},
    )
    apply_wes_panel_to_job_params(job)
    # Gene-list panels with narrow_with_panel_bed=false do not merge disease_bed into params.
    assert job.params.get("panel_interpretation_genes")
    assert job.params.get("panel_filter_after_analysis") is True
    assert job.params.get("disease_bed") is None


def test_api_payload_shape():
    rows = panels_for_api_response()
    assert isinstance(rows, list)
    if rows:
        r = rows[0]
        assert "id" in r and "label" in r and "category" in r


def test_interpretation_gene_set_merges_extra():
    job = Job(
        order_id="x",
        service_code="carrier_screening",
        sample_name="x",
        work_dir="01",
        params={
            "panel_interpretation_genes": ["BRCA1"],
            "interpretation_genes_extra": "TP53, brca2",
            "panel_filter_after_analysis": True,
        },
    )
    s = interpretation_gene_set_for_job(job)
    assert s == {"BRCA1", "TP53", "BRCA2"}


def test_should_apply_interpretation_post_filter_respects_false():
    job = Job(
        order_id="x",
        service_code="carrier_screening",
        sample_name="x",
        work_dir="01",
        params={
            "panel_filter_after_analysis": False,
            "panel_interpretation_genes": ["BRCA1"],
        },
    )
    assert should_apply_interpretation_post_filter(job) is False


def test_split_genes_trailing_period_and_multiline():
    text = (
        "ABCA3, ABCC8, CFTR,\n"
        "SMN1, ZNF469."
    )
    g = split_genes_from_text(text)
    assert "ZNF469" in g
    assert "ZNF469." not in "".join(g)
    assert g[-1] == "ZNF469"
    assert len(g) == 5


def test_filter_variants_by_interpretation_genes():
    ann = [{"gene": "CFTR"}, {"gene": "BRCA1"}]
    acmg = [{"a": 1}, {"a": 2}]
    a2, c2 = _filter_variants_by_interpretation_genes(ann, acmg, {"BRCA1"})
    assert len(a2) == 1 and a2[0]["gene"] == "BRCA1" and c2 == [{"a": 2}]
