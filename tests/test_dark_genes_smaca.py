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
