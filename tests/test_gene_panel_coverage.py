"""Gene vs panel coverage report (portal + GET /order/{id}/gene-coverage/{gene})."""

import pytest

from app.models import Job
from app.services.gene_panel_coverage import build_gene_panel_coverage_report, overlap_bp
from app.services.wes_panels import apply_wes_panel_to_job_params


def test_overlap_bp_math():
    assert overlap_bp(100, 200, 150, 180) == 30
    assert overlap_bp(100, 200, 0, 50) == 0
    assert overlap_bp(0, 1000, 100, 200) == 100


def test_invalid_gene_symbol():
    job = Job(
        order_id="x",
        service_code="carrier_screening",
        sample_name="s",
        work_dir="00",
        params={},
    )
    with pytest.raises(ValueError):
        build_gene_panel_coverage_report(job, "BAD SYM")


def test_report_bed_and_interpretation(wes_test_catalog):
    job = Job(
        order_id="ord_cov_1",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={"wes_panel_id": "test_carrier_panel"},
    )
    apply_wes_panel_to_job_params(job)
    r = build_gene_panel_coverage_report(job, "ATAD3")
    assert r["gene"] == "ATAD3"
    assert r["wes_panel_id"] == "test_carrier_panel"
    assert r["in_interpretation_set"] is True
    assert len(r["bed_regions"]) >= 1
    assert r["bed_regions"][0]["chrom"] == "chr1"
    assert r["total_target_bp"] > 0


def test_bed_col4_twist_style_semicolon_gene(tmp_path):
    """Twist-style column 4: GENE;NM_...;ENST... — match HGNC token, not full string."""
    bed = tmp_path / "twist_like.bed"
    bed.write_text(
        "chr12\t102838398\t102840000\tPAH;NM_000277.3;ENST00000307000.7\n",
        encoding="utf-8",
    )
    job = Job(
        order_id="ord_twist_bed",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={"disease_bed": str(bed.resolve())},
    )
    r = build_gene_panel_coverage_report(job, "PAH")
    assert r["gene"] == "PAH"
    assert len(r["bed_regions"]) == 1
    assert r["bed_regions"][0]["chrom"] == "chr12"
    assert r["total_target_bp"] > 0


def test_per_gene_depth_sidecar(wes_test_catalog, tmp_path):
    qc = tmp_path / "qc"
    qc.mkdir(parents=True)
    (qc / "qc_gene_coverage.tsv").write_text("gene\tmean_coverage\nATAD3\t88.5\n", encoding="utf-8")
    job = Job(
        order_id="ord_cov_2",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={"wes_panel_id": "test_carrier_panel"},
        analysis_dir=str(tmp_path),
    )
    apply_wes_panel_to_job_params(job)
    r = build_gene_panel_coverage_report(job, "ATAD3")
    assert r["per_gene_depth"] is not None
    assert r["per_gene_depth"]["mean_coverage"] == pytest.approx(88.5)


def test_api_route_exists():
    """FastAPI route table includes gene-coverage."""
    from app.main import app

    paths = [getattr(r, "path", None) for r in app.routes]
    assert any(p and "/gene-coverage/" in p for p in paths)
