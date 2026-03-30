"""Unit tests for SMAca SNP C/T extraction (import module file directly; no pydantic)."""

import importlib.util
import os

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_dark_genes():
    path = os.path.join(_ROOT, "app", "services", "carrier_screening", "dark_genes.py")
    spec = importlib.util.spec_from_file_location("dark_genes_smaca_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


_dg = _load_dark_genes()
_smaca_extract_snp_ct_counts = _dg._smaca_extract_snp_ct_counts
_try_smaca_kv_html = _dg._try_smaca_kv_html
_parse_smaca_ct_sidecar_text = _dg._parse_smaca_ct_sidecar_text
_inject_smaca_ct_line_into_detailed_text = _dg._inject_smaca_ct_line_into_detailed_text


def test_smaca_extract_ct_counts_variants():
    assert _smaca_extract_snp_ct_counts("CT_counts=10, 20") == ("10", "20")
    assert _smaca_extract_snp_ct_counts("C_T=51,49") == ("51", "49")
    assert _smaca_extract_snp_ct_counts("SMN_C_T=3, 4") == ("3", "4")
    assert _smaca_extract_snp_ct_counts("C_reads=100\nT_reads=200") == ("100", "200")
    assert _smaca_extract_snp_ct_counts("g.27134T>G AD=12,15") == ("12", "15")
    assert _smaca_extract_snp_ct_counts("C_Ratio=0.5") is None


def test_try_smaca_kv_orders_cov_fraction_before_counts():
    html = _try_smaca_kv_html(
        "SMAca CHECK (Silent Carrier + Coverage)",
        "C_Ratio=0.480\nC_T=40,44\nSMN1_CN=3",
    )
    assert html
    assert "SMN1 cov fraction" in html
    assert "0.480" in html
    assert "40 / 44" in html
    pos_cov = html.find("SMN1 cov fraction")
    pos_ct = html.find("40 / 44")
    assert pos_cov < pos_ct


def test_try_smaca_kv_cn_est_unified_pipeline():
    """Unified pipeline emits SMN1_CN_est / SMN2_CN_est (not SMN1_CN / SMN2_CN)."""
    html = _try_smaca_kv_html(
        "SMAca CHECK (exon7/c840: Pi_b cov_SMN*_b; silent carrier: g.27134 SilentCarrier):",
        "SMN1_CN_est=2\nSMN2_CN_est=2\nC_Ratio=0.480\nCov(1,2)=2.82,3.06\nSilentCarrier=False\n",
    )
    assert html
    assert "SMN1 CNV (est.)" in html
    assert "SMN2 CNV (est.)" in html
    assert ">2<" in html.replace(" ", "")


def test_parse_smaca_ct_sidecar_text():
    assert _parse_smaca_ct_sidecar_text("C_T=40, 44\n") == "C_T=40, 44"
    assert _parse_smaca_ct_sidecar_text("# comment\nCT_counts=10, 20") == "CT_counts=10, 20"
    assert _parse_smaca_ct_sidecar_text("26 16") == "C_T=26,16"


def test_inject_smaca_ct_after_c_ratio():
    body = "SMN1_CN=3\nC_Ratio=0.48\nSMN2_CN=3\n"
    out = _inject_smaca_ct_line_into_detailed_text(body, "C_T=26,16")
    assert "C_T=26,16" in out
    assert out.index("C_Ratio") < out.index("C_T=26,16")


def test_inject_idempotent_when_c_t_present():
    body = "C_T=1,2\nC_Ratio=0.5\n"
    assert _inject_smaca_ct_line_into_detailed_text(body, "C_T=9,9") == body
