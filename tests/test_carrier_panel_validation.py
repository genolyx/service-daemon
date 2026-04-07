"""Carrier screening: interpretation panel required for standard orders (not extended programs)."""

from app.services.carrier_screening.plugin import CarrierScreeningPlugin
from app.services.whole_exome import WholeExomePlugin


def test_standard_carrier_requires_wes_panel_id():
    pl = CarrierScreeningPlugin()
    ok, msg = pl.validate_params({"carrier": {"test_category": "standard_carrier"}})
    assert ok is False
    assert "wes_panel_id" in msg.lower()


def test_extended_program_skips_panel_requirement():
    pl = CarrierScreeningPlugin()
    ok, _ = pl.validate_params(
        {"carrier": {"test_category": "other", "other_test_type": "Exome"}}
    )
    assert ok is True


def test_known_panel_passes(wes_test_catalog):
    pl = CarrierScreeningPlugin()
    ok, msg = pl.validate_params(
        {
            "carrier": {"test_category": "standard_carrier"},
            "wes_panel_id": "test_carrier_panel",
        }
    )
    assert ok is True, msg


def test_unknown_panel_fails():
    pl = CarrierScreeningPlugin()
    ok, msg = pl.validate_params(
        {
            "carrier": {"test_category": "standard_carrier"},
            "wes_panel_id": "___no_such_panel___",
        }
    )
    assert ok is False
    assert "unknown" in msg.lower()


def test_whole_exome_validate_params_accepts_strict_kwarg():
    """save_order passes strict=False; WholeExome must accept it (no TypeError → 500)."""
    pl = WholeExomePlugin()
    ok, _ = pl.validate_params({}, strict=False)
    assert ok is True
