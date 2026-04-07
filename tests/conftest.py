"""Shared pytest fixtures."""

from pathlib import Path

import pytest


@pytest.fixture
def wes_test_catalog(monkeypatch):
    """Point the bundled WES panel catalog at a one-panel fixture for tests."""
    from app.config import settings

    p = Path(__file__).resolve().parent / "fixtures" / "wes_panel_catalog_test.json"
    monkeypatch.setattr(settings, "wes_panels_json", str(p))
