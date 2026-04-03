"""
Health Screening — same gx-exome Nextflow pipeline and review stack as whole exome /
carrier screening. Primary (interpretation) panel is required on submit/run (strict validation),
same as standard carrier screening; draft save may omit the panel until it is chosen.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple

from ..carrier_screening.plugin import CarrierScreeningPlugin
from ..whole_exome import WholeExomePlugin


class HealthScreeningPlugin(WholeExomePlugin):
    @property
    def service_code(self) -> str:
        return "health_screening"

    @property
    def display_name(self) -> str:
        return "Health Screening"

    def validate_params(self, params: Dict[str, Any], strict: bool = True) -> Tuple[bool, str]:
        """Like carrier_screening: ``wes_panel_id`` required when ``strict=True`` (submit/start)."""
        return CarrierScreeningPlugin.validate_params(self, params, strict)


def create_plugin():
    return HealthScreeningPlugin()
