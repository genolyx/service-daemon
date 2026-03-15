"""
Import test for all service-daemon modules.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

print("Testing imports...")

# Core modules
from app.config import settings
print(f"  [OK] config (enabled_services: {settings.enabled_service_list})")

from app.models import (
    Job, OrderStatus, ServiceCode, OrderSubmitRequest, OrderSubmitResponse,
    OrderStatusResponse, NotificationResult, OutputFile, QueueSummary,
    ReportGenerateRequest, ReportGenerateResponse,
)
print(f"  [OK] models (OrderStatus values: {[s.value for s in OrderStatus]})")

from app.auth_client import AuthClient
print(f"  [OK] auth_client")

from app.platform_client import PlatformClient
print(f"  [OK] platform_client")

from app.queue_manager import QueueManager, get_queue_manager
print(f"  [OK] queue_manager")

from app.runner import PipelineRunner
print(f"  [OK] runner")

# Service plugins
from app.services.base import ServicePlugin
print(f"  [OK] services.base")

from app.services import load_plugins, get_plugin, list_service_codes
print(f"  [OK] services registry")

# Load plugins
load_plugins(["carrier_screening", "nipt", "sgnipt"])
codes = list_service_codes()
print(f"  [OK] Loaded plugins: {codes}")

# Verify each plugin
for code in codes:
    plugin = get_plugin(code)
    assert plugin is not None, f"Plugin not found: {code}"
    assert plugin.service_code == code, f"Service code mismatch: {plugin.service_code} != {code}"
    print(f"  [OK] Plugin '{code}' -> {plugin.display_name}")
    print(f"       Progress stages: {plugin.get_progress_stages()}")

# Test carrier_screening specific components
from app.services.carrier_screening import CarrierScreeningPlugin
print(f"  [OK] carrier_screening.CarrierScreeningPlugin")

from app.services.carrier_screening.annotator import (
    ClinVarAnnotator, GnomADAnnotator, ClinGenAnnotator,
    MANEAnnotator, GeneIntervalAnnotator,
    HPOAnnotator, CuratedVariantDB, HGMDAnnotator, DiseaseGeneMapper,
    VariantAnnotator,
)
print(f"  [OK] carrier_screening.annotator (all classes)")

from app.services.carrier_screening.vcf_parser import (
    clean_vcf_remove_formats, get_annotation_layout,
    extract_variant_info, get_sample_metrics,
    load_bed_regions, variant_in_bed,
    VariantFilterConfig, apply_clinvar_filter,
    parse_vcf_variants,
)
print(f"  [OK] carrier_screening.vcf_parser (all functions)")

from app.services.carrier_screening.review import (
    extract_qc_summary, generate_result_json, generate_variants_tsv,
)
print(f"  [OK] carrier_screening.review (all functions)")

from app.services.carrier_screening.report import (
    generate_report_json, generate_report_pdf,
)
print(f"  [OK] carrier_screening.report (all functions)")

from app.services.carrier_screening.acmg import classify_variant
print(f"  [OK] carrier_screening.acmg")

# Test Job creation
job = Job(
    order_id="TEST-001",
    service_code="carrier_screening",
    sample_name="SAMPLE_001",
    work_dir="260101",
    fastq_dir="/data/fastq/260101/SAMPLE_001",
    analysis_dir="/data/analysis/carrier_screening/260101/SAMPLE_001",
    output_dir="/data/output/carrier_screening/260101/SAMPLE_001",
    log_dir="/data/log/carrier_screening/260101/SAMPLE_001"
)
print(f"  [OK] Job created: {job.order_id} ({job.service_code})")

# Test ReportGenerateRequest
report_req = ReportGenerateRequest(
    confirmed_variants=[{"variant_id": "VAR_0001", "gene": "CFTR"}],
    reviewer_info={"name": "Dr. Test", "id": "reviewer_001"},
    patient_info={"name": "John Doe", "dob": "1990-01-01"},
    partner_info={"name": "Jane Doe", "dob": "1991-02-02"},
    languages=["EN", "CN"],
)
print(f"  [OK] ReportGenerateRequest created: {len(report_req.confirmed_variants)} variants, languages={report_req.languages}")

# Test QueueManager
import asyncio

async def test_queue():
    qm = QueueManager(max_concurrent=2)
    pos = await qm.enqueue(job)
    print(f"  [OK] QueueManager: enqueued at position {pos}")
    summary = qm.get_summary()
    print(f"  [OK] QueueSummary: queued={summary.total_queued}, running={summary.total_running}")

asyncio.run(test_queue())

# Test FastAPI app import
from app.main import app
print(f"  [OK] FastAPI app: {app.title} v{app.version}")
routes = [r.path for r in app.routes if hasattr(r, 'path')]
print(f"  [OK] Routes: {routes}")

print("\n✅ All imports and basic tests passed!")
