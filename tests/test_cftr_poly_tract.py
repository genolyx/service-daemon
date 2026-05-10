"""CFTR poly-T / poly-TG parsing from dark-gene report blobs."""

import importlib.util
from pathlib import Path

_DG_MOD = None


def _dg():
    global _DG_MOD
    if _DG_MOD is None:
        p = Path(__file__).resolve().parent.parent / "app/services/carrier_screening/dark_genes.py"
        spec = importlib.util.spec_from_file_location("_dg_cftr", p)
        _DG_MOD = importlib.util.module_from_spec(spec)
        assert spec.loader
        spec.loader.exec_module(_DG_MOD)
    return _DG_MOD


def extract_cftr_poly_tract_calls(text: str):
    return _dg().extract_cftr_poly_tract_calls(text)


def test_extract_kv_poly_t_and_tg():
    raw = """
    SAMPLE: test
    poly_t=9T/9T
    poly_tg=TG11/TG12
    """
    d = extract_cftr_poly_tract_calls(raw)
    assert d.get("T") == "9T/9T"
    assert d.get("TG") == "TG11/TG12"


def test_extract_labeled_lines():
    raw = """
    CFTR poly-T: 5T/7T
    poly-TG: TG11/TG12
    """
    d = extract_cftr_poly_tract_calls(raw)
    assert d.get("T") == "5T/7T"
    assert d.get("TG") == "TG11/TG12"


def test_extract_tg_line():
    raw = "TG tract:\tTG11/TG11\n"
    d = extract_cftr_poly_tract_calls(raw)
    assert d.get("TG") == "TG11/TG11"


def test_extract_ivs_headers():
    raw = """
    IVS8 poly-T: 9T/7T
    IVS9 poly-TG: (TG)11/(TG)12
    """
    d = extract_cftr_poly_tract_calls(raw)
    assert "T" in d
    assert "TG" in d


def test_extract_from_dark_genes_sections_only():
    """CFTR block may live only in structured sections (same as Dark genes tab)."""
    extract_fn = _dg().extract_cftr_poly_tract_calls_from_dark_genes

    blk = {
        "summary_text": "",
        "detailed_text": "",
        "detailed_sections": [
            {
                "title": "CFTR screening",
                "body": "Poly T:\t5T/7T\nPoly TG:\tTG10/TG11\n",
            }
        ],
    }
    d = extract_fn(blk)
    assert d.get("T") == "5T/7T"
    assert d.get("TG") == "TG10/TG11"


def test_spaced_poly_t_label():
    d = extract_cftr_poly_tract_calls("Poly T: 9T/9T\n")
    assert d.get("T") == "9T/9T"


def test_cftr_context_guess():
    d = extract_cftr_poly_tract_calls(
        "CFTR\nIVS8 genotype  5T/7T  (informative)\n"
    )
    assert d.get("T") == "5T/7T"


def test_expansion_hunter_ivs9_repcn():
    mod = _dg()
    sample = """
  Raw EH REPCN  : CFTR_TG=11/11  CFTR_polyT=7/7
  Per-allele    : (TG)11-7T [BENIGN] | (TG)11-7T [BENIGN]
"""
    eh = mod.parse_cftr_expansion_hunter_ivs9(sample)
    assert eh is not None
    assert eh["raw_poly_t"] == "7/7"
    assert eh["raw_tg"] == "11/11"
    assert eh["risk_level"] == "low"
    d = mod.extract_cftr_poly_tract_calls(sample)
    assert d.get("T") == "7T/7T"
    assert d.get("TG") == "TG11/TG11"


def test_expansion_hunter_flags_5t_tg12():
    mod = _dg()
    eh = mod.parse_cftr_expansion_hunter_ivs9("CFTR_polyT=5/7 CFTR_TG=11/12\n")
    assert eh["risk_level"] == "high"
    assert any("5T" in r for r in eh["risk_reasons"])
    assert any("12" in r for r in eh["risk_reasons"])


def test_expansion_hunter_alternate_spaced_labels():
    mod = _dg()
    eh = mod.parse_cftr_expansion_hunter_ivs9("CFTR polyT=7/7 CFTR TG=11/11\n")
    assert eh is not None
    assert eh["raw_poly_t"] == "7/7"
    assert eh["raw_tg"] == "11/11"
    assert eh["risk_level"] == "low"


def test_expansion_hunter_flags_tg13():
    mod = _dg()
    eh = mod.parse_cftr_expansion_hunter_ivs9("CFTR_polyT=7/7 CFTR_TG=11/13\n")
    assert eh["risk_level"] == "high"
    assert any("13" in r and "TG" in r for r in eh["risk_reasons"])


def test_pdf_cftr_eh_section_emits_only_poly_t_and_tg_rows():
    mod = _dg()
    body = """Raw EH REPCN  : CFTR_TG=11/11  CFTR_polyT=7/7
Per-allele    : foo
Risk assessment : bar
"""
    html_out = mod._section_body_portal_html_for_pdf({"title": "CFTR", "body": body, "kind": "normal"})
    assert "Poly-T" in html_out
    assert "(EH REPCN)" not in html_out
    assert "7 / 7" in html_out
    assert "11 / 11" in html_out
    assert "Per-allele" not in html_out
    assert "Raw EH REPCN" not in html_out
    assert "Risk" not in html_out


def test_pdf_cftr_benign_eh_omitted_from_customer_pdf_even_if_approved():
    """EH 7/7 + TG11/11 must not land on the customer PDF when explicitly approved."""
    mod = _dg()
    pdf = mod.dark_genes_for_pdf(
        {
            "status": "found",
            "detailed_sections": [
                {
                    "title": "CFTR IVS8 (Expansion Hunter)",
                    "body": "Raw EH REPCN  : CFTR_TG=11/11  CFTR_polyT=7/7\n",
                    "kind": "normal",
                },
            ],
            "section_reviews": [{"approved": True, "risk": "low", "notes": ""}],
        }
    )
    assert not pdf


def test_dark_genes_pdf_supplemental_high_risk_row():
    mod = _dg()
    pdf = mod.dark_genes_for_pdf(
        {
            "status": "found",
            "detailed_sections": [
                {
                    "title": "CFTR adjunct",
                    "body": "CFTR_polyT=5/7 CFTR_TG=11/11\n",
                    "kind": "normal",
                },
            ],
            "section_reviews": [{"approved": True, "risk": "high", "notes": "Lab confirm"}],
        }
    )
    assert pdf.get("report_detailed_html")
    assert pdf.get("supplemental_review_high_risk") is True
    rows = pdf.get("supplemental_summary_findings") or []
    assert len(rows) == 1
    assert rows[0].get("finding_source") == "dark_genes_supplemental"
    assert rows[0].get("gene") == "CFTR"
    assert "High priority" in (rows[0].get("classification") or "")


def test_cftr_call_negative_filtered():
    assert not _dg().cftr_tract_call_is_actionable("negative")
    assert not _dg().cftr_tract_call_is_actionable("N/A")
    assert _dg().cftr_tract_call_is_actionable("9T/9T")


def test_section_without_cftr_title_still_scanned():
    mod = _dg()
    blk = {
        "summary_text": "",
        "detailed_text": "",
        "detailed_sections": [
            {
                "title": "Adjunct intronic repeats",
                "body": "poly_t=9T/7T\n",
            },
        ],
    }
    d = mod.extract_cftr_poly_tract_calls_from_dark_genes(blk)
    assert d.get("T") == "9T/7T"


def test_pdf_detailed_sections_core_always_included_low_risk_adjunct_not_auto_approved():
    """Core blocks always on PDF without approval; low-risk non-core still needs explicit approval."""
    mod = _dg()
    sections = [
        {"title": "Spinal Muscular Atrophy", "body": "Carrier screen detail.\n", "kind": "normal"},
        {
            "title": "HBA Analysis (Alpha Thalassemia - Dosage)",
            "body": "Dosage detail.\n",
            "kind": "normal",
        },
        {"title": "CFTR screening", "body": "Poly T:\t9T/9T\n", "kind": "normal"},
    ]
    reviews = [
        {"approved": False, "notes": "", "risk": "low"},
        {"approved": False, "notes": "", "risk": "low"},
        {"approved": False, "notes": "", "risk": "low"},
    ]
    html_out = mod.detailed_sections_to_pdf_html(sections, reviews, filter_by_approval=True)
    assert "Spinal Muscular Atrophy" in html_out
    assert "Alpha Thalassemia" in html_out
    assert "CFTR" not in html_out


def test_pdf_detailed_sections_high_risk_adjunct_excluded_when_not_effectively_approved():
    """High-risk non-core stays off the customer PDF until explicitly approved + stored high."""
    mod = _dg()
    sections = [
        {"title": "Spinal Muscular Atrophy", "body": "Carrier screen detail.\n", "kind": "normal"},
        {
            "title": "Large SVs (Manta/gCNV - rest)",
            "body": "WARNING: possible artifact\n",
            "kind": "normal",
        },
    ]
    reviews = [
        {"approved": False, "notes": "", "risk": "low"},
        {"approved": False, "notes": "", "risk": "low"},
    ]
    html_out = mod.detailed_sections_to_pdf_html(sections, reviews, filter_by_approval=True)
    assert "Spinal Muscular Atrophy" in html_out
    assert "Structural" not in html_out and "CNV" not in html_out


def test_pdf_dark_genes_high_risk_core_omitted_until_approved():
    """High-risk core block does not appear on PDF (and no supplemental row) until approved."""
    mod = _dg()
    pdf = mod.dark_genes_for_pdf(
        {
            "status": "found",
            "detailed_sections": [
                {"title": "Spinal Muscular Atrophy", "body": "Lab notes.\n", "kind": "normal"},
            ],
            "section_reviews": [{"approved": False, "risk": "high", "notes": "review pending"}],
        }
    )
    assert not pdf


def test_pdf_cah_high_unapproved_omitted_from_customer_html():
    mod = _dg()
    sections = [
        {
            "title": "CYP21A2 analysis (CAH - dosage)",
            "body": "Possible deletion.\n",
            "kind": "normal",
        },
    ]
    reviews = [{"approved": False, "risk": "high", "notes": ""}]
    html_out = mod.detailed_sections_to_pdf_html(sections, reviews, filter_by_approval=True)
    assert not html_out.strip()
