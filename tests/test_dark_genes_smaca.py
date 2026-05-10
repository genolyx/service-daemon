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


_infer_high = _dg._infer_pipeline_section_high_risk
_align_rev = _dg.align_section_reviews


def test_cah_possible_deletion_flags_pipeline_high_risk():
    sec = {
        "title": "CYP21A2 analysis (CAH - dosage)",
        "body": "Dosage supports a possible deletion in the CYP21A2 / paralog region.\n",
        "kind": "normal",
    }
    assert _infer_high(sec) is True


def test_cah_negated_possible_deletion_not_flagged():
    sec = {
        "title": "CYP21A2 analysis (CAH - dosage)",
        "body": "This pattern is not a possible deletion after re-check.\n",
        "kind": "normal",
    }
    assert _infer_high(sec) is False


def test_cah_paralog_deletion_kv_possible_flags_high_risk():
    sec = {
        "title": "CYP21A2 analysis (CAH - dosage)",
        "body": "paralog_deletion=possible\n",
        "kind": "normal",
    }
    assert _infer_high(sec) is True


def test_possible_deletion_non_cah_section_not_flagged():
    sec = {"title": "Large SVs (other)", "body": "possible deletion noted elsewhere\n", "kind": "normal"}
    assert _infer_high(sec) is False


def test_align_section_reviews_default_risk_high_for_cah_possible_deletion():
    sections = [
        {
            "title": "CYP21A2 analysis (CAH - dosage)",
            "body": "Summary: possible deletion\n",
            "kind": "normal",
        },
    ]
    out = _align_rev(None, 1, sections)
    assert out[0]["risk"] == "high"
    assert out[0]["approved"] is False


def test_align_core_section_cftr_lines_in_body_still_low_auto_approved():
    """CFTR EH lines inside an HBA core block must not disable core low-risk auto-approve."""
    sections = [
        {
            "title": "HBA Analysis (Alpha Thalassemia - Dosage)",
            "body": "CFTR_polyT=7/7 CFTR_TG=11/11\nhba1=2\n",
            "kind": "normal",
        },
    ]
    out = _align_rev(None, 1, sections)
    assert out[0]["risk"] == "low"
    assert out[0]["approved"] is True


def test_cftr_negative_title_never_low_risk_auto_effective_approved():
    sec = {
        "title": "CFTR IVS9 adjunct",
        "body": "CFTR_polyT=7/7 CFTR_TG=11/11\n",
        "kind": "normal",
    }
    rev = {"approved": False, "risk": "low", "notes": ""}
    assert _dg.cftr_supplementary_section_excludes_low_risk_auto_approve(sec) is True
    assert _dg.effective_approved_for_dark_genes_section(rev, sec) is False


def test_cftr_supplementary_detects_poly_t_body_without_cftr_title():
    sec = {
        "title": "Adjunct intronic repeats",
        "body": "poly_t=9T/9T\n poly_tg=TG11/TG11\n",
        "kind": "normal",
    }
    assert _dg.cftr_supplementary_section_excludes_low_risk_auto_approve(sec) is True


def test_align_section_reviews_low_risk_auto_approved_core_only():
    sections = [
        {
            "title": "CYP21A2 analysis (CAH - dosage)",
            "body": "Normal copy-number result, no warnings.\n",
            "kind": "normal",
        },
    ]
    out = _align_rev(None, 1, sections)
    assert out[0]["risk"] == "low"
    assert out[0]["approved"] is True


def test_align_section_reviews_low_risk_non_core_not_auto_approved():
    sections = [
        {
            "title": "CFTR screening",
            "body": "Poly T:\t9T/9T\n",
            "kind": "normal",
        },
    ]
    out = _align_rev(None, 1, sections)
    assert out[0]["risk"] == "low"
    assert out[0]["approved"] is False


def test_effective_approved_true_when_effective_low_core_even_if_disk_flag_false():
    sec = {
        "title": "CYP21A2 analysis (CAH - dosage)",
        "body": "Benign.\n",
        "kind": "normal",
    }
    rev = {"approved": False, "risk": "low", "notes": ""}
    assert _dg.effective_approved_for_dark_genes_section(rev, sec) is True


def test_effective_approved_cftr_molecular_benign_never_customer_pdf_eligible():
    """Benign CFTR tracts are not customer-PDF-eligible even when the portal marks approved."""
    sec = {"title": "CFTR screening", "body": "Poly T:\t9T/9T\n", "kind": "normal"}
    assert _dg.effective_approved_for_dark_genes_section(
        {"approved": False, "risk": "low", "notes": ""}, sec
    ) is False
    assert _dg.effective_approved_for_dark_genes_section(
        {"approved": True, "risk": "low", "notes": ""}, sec
    ) is False


def test_effective_approved_cftr_molecular_high_non_core_requires_portal_approved():
    sec = {
        "title": "CFTR screening",
        "body": "CFTR_polyT=5/7 CFTR_TG=11/11\n",
        "kind": "normal",
    }
    assert _dg.effective_approved_for_dark_genes_section(
        {"approved": False, "risk": "low", "notes": ""}, sec
    ) is False
    assert _dg.effective_approved_for_dark_genes_section(
        {"approved": True, "risk": "low", "notes": ""}, sec
    ) is True


def test_effective_risk_overrides_stale_disk_low_when_cah_prose_high():
    sec = {
        "title": "CYP21A2 analysis (CAH - dosage)",
        "body": "Call: possible deletion on paralog background.\n",
        "kind": "normal",
    }
    rev = {"approved": True, "notes": "", "risk": "low"}
    assert _dg.effective_risk_for_section(rev, sec) == "high"


def test_align_upgrades_stored_low_when_inference_high_and_clears_approval():
    sections = [
        {
            "title": "CYP21A2 analysis (CAH - dosage)",
            "body": "Interpretation suggests possible deletion.\n",
            "kind": "normal",
        },
    ]
    prev = [{"approved": True, "notes": "", "risk": "low"}]
    out = _align_rev(prev, 1, sections)
    assert out[0]["risk"] == "high"
    assert out[0]["approved"] is False


def test_align_clears_approval_when_inference_high_and_disk_risk_missing():
    sections = [
        {
            "title": "CYP21A2 analysis (CAH - dosage)",
            "body": "Possible deletion.\n",
            "kind": "normal",
        },
    ]
    prev = [{"approved": True, "notes": "", "risk": None}]
    out = _align_rev(prev, 1, sections)
    assert out[0]["risk"] == "high"
    assert out[0]["approved"] is False


def test_effective_approved_false_when_inferred_high_but_explicit_risk_not_high():
    sec = {
        "title": "CYP21A2 analysis (CAH - dosage)",
        "body": "Possible deletion.\n",
        "kind": "normal",
    }
    rev_lo = {"approved": True, "risk": "low", "notes": ""}
    rev_miss = {"approved": True, "notes": ""}
    assert _dg.effective_approved_for_dark_genes_section(rev_lo, sec) is False
    assert _dg.effective_approved_for_dark_genes_section(rev_miss, sec) is False
    rev_hi = {"approved": True, "risk": "high", "notes": ""}
    assert _dg.effective_approved_for_dark_genes_section(rev_hi, sec) is True


def test_supplemental_high_risk_skips_stale_approve_without_explicit_high_risk():
    pdf = _dg.dark_genes_for_pdf(
        {
            "status": "found",
            "detailed_sections": [
                {
                    "title": "CYP21A2 analysis (CAH - dosage)",
                    "body": "Possible deletion.\n",
                    "kind": "normal",
                },
            ],
            "section_reviews": [{"approved": True, "risk": "low", "notes": ""}],
        }
    )
    assert not pdf.get("supplemental_summary_findings")


def test_cah_title_bare_cah_dosage_matches_for_deletion_prose():
    sec = {
        "title": "CAH — Dosage analysis",
        "body": "Possible deletion.\n",
        "kind": "normal",
    }
    assert _infer_high(sec) is True
