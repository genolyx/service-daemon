"""
Whole Exome service — same pipeline and review stack as carrier screening, but
``wes_panel_id`` (Primary / interpretation) is optional: full exome can be reported
when no panel is selected.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple

from ..carrier_screening.plugin import CarrierScreeningPlugin, _resolve_carrier_wes_panel_id


class WholeExomePlugin(CarrierScreeningPlugin):
    @property
    def service_code(self) -> str:
        return "whole_exome"

    @property
    def display_name(self) -> str:
        return "Whole Exome"

    def validate_params(self, params: Dict[str, Any]) -> Tuple[bool, str]:
        carrier = (params or {}).get("carrier") or {}
        if carrier.get("reuse_prior_pipeline_outputs"):
            pid = (carrier.get("prior_order_id") or "").strip()
            if not pid:
                return (
                    False,
                    "prior_order_id is required when reuse_prior_pipeline_outputs is true",
                )
        wid = _resolve_carrier_wes_panel_id(params)
        if wid:
            from ..wes_panels import get_panel_by_id, resolve_panel_interpretation_genes

            panel = get_panel_by_id(wid)
            if not panel:
                return (
                    False,
                    f"Unknown wes_panel_id={wid!r} — add it to the WES panel catalog or fix the id.",
                )
            genes = resolve_panel_interpretation_genes(panel)
            if not genes:
                return (
                    False,
                    f"Panel {wid!r} resolves to zero interpretation genes "
                    "(needs interpretation_genes, interpretation_genes_file, or "
                    "interpretation_genes_from_disease_bed + disease_bed).",
                )
        return True, ""


def create_plugin():
    return WholeExomePlugin()
