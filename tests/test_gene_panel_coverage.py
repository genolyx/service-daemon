"""Gene vs panel coverage report (portal + GET /order/{id}/gene-coverage/{gene})."""

import pytest

from app.config import settings
from app.models import Job
from app.services import gene_panel_coverage as gpc
from app.services.gene_panel_coverage import (
    _exon_coverage_rows,
    _exons_for_gene_from_gff,
    _merge_bed_regions_half_open,
    build_gene_panel_coverage_report,
    overlap_bp,
)
from app.services.wes_panels import apply_wes_panel_to_job_params


@pytest.fixture
def twist_bed_pa(tmp_path):
    p = tmp_path / "targets.bed"
    p.write_text(
        "chr12\t102838398\t102840000\tPAH;NM_000277.3;ENST00000307000.7\n",
        encoding="utf-8",
    )
    return p


def test_merge_bed_regions_adjacent_and_gap():
    raw = [
        {"chrom": "chr1", "start": 10, "end": 11, "name": "G", "length_bp": 1},
        {"chrom": "chr1", "start": 11, "end": 20, "name": "G", "length_bp": 9},
        {"chrom": "chr1", "start": 100, "end": 105, "name": "G", "length_bp": 5},
    ]
    m = _merge_bed_regions_half_open(raw)
    assert len(m) == 2
    assert m[0]["start"] == 10 and m[0]["end"] == 20 and m[0]["length_bp"] == 10
    assert m[1]["start"] == 100 and m[1]["end"] == 105


def test_report_merges_multiple_twist_rows(tmp_path, monkeypatch):
    tb = tmp_path / "t.bed"
    tb.write_text(
        "chr1\t10\t11\tFOO\nchr1\t11\t20\tFOO\nchr1\t100\t105\tFOO\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(gpc, "_twist_exome_targets_bed", lambda _j: str(tb.resolve()))
    job = Job(
        order_id="ord_merge",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={},
    )
    r = build_gene_panel_coverage_report(job, "FOO")
    assert r["twist_raw_interval_count"] == 3
    assert len(r["bed_regions"]) == 2
    assert r["total_target_bp"] == 10 + 5


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


def test_report_uses_twist_bed_only_ignore_carrier_beds(tmp_path, monkeypatch, wes_test_catalog):
    """Interpretation panel is genes only; depth intervals come from Twist exome BED, not disease/backbone params."""
    tb = tmp_path / "twist.bed"
    tb.write_text("chr1\t100\t500\tATAD3;NM_001\n", encoding="utf-8")
    monkeypatch.setattr(gpc, "_twist_exome_targets_bed", lambda _j: str(tb.resolve()))

    job = Job(
        order_id="ord_cov_1",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={
            "wes_panel_id": "test_carrier_panel",
            "disease_bed": "/dev/null/should_not_be_used.bed",
            "backbone_bed": "/dev/null/also_ignored.bed",
        },
    )
    apply_wes_panel_to_job_params(job)
    r = build_gene_panel_coverage_report(job, "ATAD3")
    assert r["gene"] == "ATAD3"
    assert r["wes_panel_id"] == "test_carrier_panel"
    assert r["in_interpretation_set"] is True
    assert r["disease_bed_path"] == str(tb.resolve())
    assert r["disease_bed_source"] == "twist_exome2.targets_bed"
    assert r["clinical_disease_bed_path"] is None
    assert r["intervals_clipped_to_capture"] is False
    assert len(r["bed_regions"]) >= 1
    assert r["bed_regions"][0]["chrom"] == "chr1"
    assert r["total_target_bp"] > 0


def test_twist_col4_semicolon_token(monkeypatch, twist_bed_pa):
    monkeypatch.setattr(gpc, "_twist_exome_targets_bed", lambda _j: str(twist_bed_pa.resolve()))
    job = Job(
        order_id="ord_twist",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={},
    )
    r = build_gene_panel_coverage_report(job, "PAH")
    assert r["gene"] == "PAH"
    assert len(r["bed_regions"]) == 1
    assert r["bed_regions"][0]["chrom"] == "chr12"
    assert r["total_target_bp"] > 0


def test_twist_exome_targets_bed_ignores_capture_panel_id(tmp_path, monkeypatch):
    """Always resolve twist-exome2 layout, not job.params.carrier.capture_panel_id."""
    bed_dir = tmp_path / "data" / "bed" / "twist-exome2"
    bed_dir.mkdir(parents=True)
    targets = bed_dir / "targets.bed"
    targets.write_text("chr12\t1\t100\tPAH\n", encoding="utf-8")

    class FakePl:
        def _resolve_capture_panel_bed(self, capture_panel_id: str):
            p = tmp_path / "data" / "bed" / capture_panel_id / "targets.bed"
            return str(p) if p.is_file() else None

        def _run_analysis_data_dir(self):
            return str(tmp_path)

    monkeypatch.setattr("app.services.get_plugin", lambda _code: FakePl())

    job = Job(
        order_id="ord_seq",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={"carrier": {"capture_panel_id": "some-other-panel"}},
    )
    assert gpc._twist_exome_targets_bed(job) == str(targets.resolve())


def test_per_gene_depth_sidecar(wes_test_catalog, tmp_path, monkeypatch):
    tb = tmp_path / "t.bed"
    tb.write_text("chr1\t1\t999\tATAD3\n", encoding="utf-8")
    monkeypatch.setattr(gpc, "_twist_exome_targets_bed", lambda _j: str(tb.resolve()))

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


def test_exons_from_gff_picks_most_exon_rich_transcript(tmp_path):
    gff = tmp_path / "m.gff3"
    gff.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\t.\texon\t101\t200\t.\t.\t.\tParent=txB;gene=FOO;exon_number=1",
                "chr1\t.\texon\t301\t400\t.\t.\t.\tParent=txA;gene=FOO;exon_number=1",
                "chr1\t.\texon\t401\t500\t.\t.\t.\tParent=txA;gene=FOO;exon_number=2",
            ]
        ),
        encoding="utf-8",
    )
    ex, tid = _exons_for_gene_from_gff(str(gff.resolve()), "FOO")
    assert tid == "txA"
    assert len(ex) == 2
    assert ex[0]["start"] == 300 and ex[0]["end"] == 400 and ex[0]["length_bp"] == 100
    assert ex[1]["start"] == 400 and ex[1]["end"] == 500


def test_exon_coverage_rows_without_mosdepth(tmp_path):
    gff = tmp_path / "m.gff3"
    gff.write_text(
        "chr1\t.\texon\t2\t3\t.\t.\t.\tParent=t1;gene_name=BAR;exon_number=1\n",
        encoding="utf-8",
    )
    ex, _tid = _exons_for_gene_from_gff(str(gff.resolve()), "BAR")
    rows = _exon_coverage_rows(ex, None)
    assert len(rows) == 1
    assert rows[0]["coverage_quality"] == "unknown"
    assert rows[0]["pct_bases_ge_20x"] is None
    assert rows[0]["start"] == 1 and rows[0]["end"] == 3 and rows[0]["length_bp"] == 2


def test_report_includes_exon_coverage_when_mane_gff_set(tmp_path, monkeypatch, wes_test_catalog):
    tb = tmp_path / "t.bed"
    tb.write_text("chr1\t1\t999\tATAD3\n", encoding="utf-8")
    gff = tmp_path / "mane.gff3"
    gff.write_text(
        "chr1\t.\texon\t10\t19\t.\t.\t.\tParent=NM_fake;gene_name=ATAD3;exon_number=1\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(gpc, "_twist_exome_targets_bed", lambda _j: str(tb.resolve()))
    monkeypatch.setattr(settings, "mane_gff", str(gff.resolve()))

    job = Job(
        order_id="ord_ex",
        service_code="carrier_screening",
        sample_name="sam",
        work_dir="00",
        params={"wes_panel_id": "test_carrier_panel"},
    )
    apply_wes_panel_to_job_params(job)
    r = build_gene_panel_coverage_report(job, "ATAD3")
    assert r.get("exon_coverage") is not None
    ec = r["exon_coverage"]
    assert ec["transcript_id"] == "NM_fake"
    assert len(ec["exons"]) >= 1
    assert ec["exons"][0]["chrom"] == "chr1"
