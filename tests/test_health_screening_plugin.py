"""Health Screening: same gx-exome stack as whole exome; primary panel required on submit."""

from app.services.health_screening import HealthScreeningPlugin


def test_health_screening_validate_params_accepts_strict_kwarg():
    pl = HealthScreeningPlugin()
    ok, _ = pl.validate_params({}, strict=False)
    assert ok is True


def test_health_screening_strict_requires_wes_panel():
    pl = HealthScreeningPlugin()
    ok, msg = pl.validate_params({"carrier": {"test_category": "standard_carrier"}}, strict=True)
    assert ok is False
    assert "wes_panel_id" in msg.lower()


def test_health_screening_service_code():
    pl = HealthScreeningPlugin()
    assert pl.service_code == "health_screening"
    assert "Health" in pl.display_name
