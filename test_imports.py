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
    OrderStatusResponse, NotificationResult, OutputFile, QueueSummary
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
from app.services.carrier_screening import (
    VCFAnnotator, ReviewPageGenerator, parse_vcf_variants,
    CarrierScreeningPlugin
)
print(f"  [OK] carrier_screening components")

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
