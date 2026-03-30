"""
API endpoint test using FastAPI TestClient.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

# Pre-load plugins (lifespan runs async, TestClient may not fully execute it)
from app.services import load_plugins
load_plugins(["carrier_screening", "nipt", "sgnipt"])

print("Testing API endpoints...\n")

# 1. Health check
resp = client.get("/health")
assert resp.status_code == 200
data = resp.json()
assert data["status"] == "healthy"
assert "carrier_screening" in data["registered_services"]
print(f"  [OK] GET /health -> {data['status']} (services: {data['registered_services']})")

# 2. Dashboard
resp = client.get("/")
assert resp.status_code == 200
data = resp.json()
assert data["service"] == "service-daemon"
print(f"  [OK] GET / -> {data['service']} v{data['version']}")

# 3. List services
resp = client.get("/services")
assert resp.status_code == 200
data = resp.json()
assert data["total"] == 3
service_codes = [s["service_code"] for s in data["services"]]
assert "carrier_screening" in service_codes
assert "nipt" in service_codes
assert "sgnipt" in service_codes
print(f"  [OK] GET /services -> {data['total']} services: {service_codes}")

# 4. Queue summary
resp = client.get("/queue/summary")
assert resp.status_code == 200
data = resp.json()
print(f"  [OK] GET /queue/summary -> queued={data['total_queued']}, running={data['total_running']}")

# 5. Submit carrier_screening order
order_payload = {
    "order_id": "CS-2026-001",
    "service_code": "carrier_screening",
    "work_dir": "260227",
    "fastq_r1_path": "/data/fastq/260227/CS-2026-001/CS-2026-001_R1.fastq.gz",
    "fastq_r2_path": "/data/fastq/260227/CS-2026-001/CS-2026-001_R2.fastq.gz",
    "params": {},
    "priority": "normal"
}
resp = client.post("/order/carrier_screening/submit", json=order_payload)
assert resp.status_code == 200
data = resp.json()
assert data["status"] == "accepted"
assert data["service_code"] == "carrier_screening"
print(f"  [OK] POST /order/carrier_screening/submit -> {data['status']} (queue_pos: {data['queue_position']})")

# 6. Submit NIPT order
nipt_payload = {
    "order_id": "NIPT-2026-001",
    "service_code": "nipt",
    "work_dir": "260227",
    "fastq_r1_path": "/data/fastq/260227/NIPT-2026-001/NIPT-2026-001_R1.fastq.gz",
    "fastq_r2_path": "/data/fastq/260227/NIPT-2026-001/NIPT-2026-001_R2.fastq.gz"
}
resp = client.post("/order/nipt/submit", json=nipt_payload)
assert resp.status_code == 200
data = resp.json()
assert data["service_code"] == "nipt"
print(f"  [OK] POST /order/nipt/submit -> {data['status']} (queue_pos: {data['queue_position']})")

# 7. Check order status
resp = client.get("/order/CS-2026-001/status")
assert resp.status_code == 200
data = resp.json()
assert data["service_code"] == "carrier_screening"
print(f"  [OK] GET /order/CS-2026-001/status -> {data['status']} ({data['progress']}%)")

# 8. Queue status with service filter
resp = client.get("/queue/status?service_code=carrier_screening")
assert resp.status_code == 200
data = resp.json()
print(f"  [OK] GET /queue/status?service_code=carrier_screening -> slots: {data['available_slots']}/{data['max_concurrent']}")

# 9. Unknown service should fail
resp = client.post("/order/unknown_service/submit", json={
    "order_id": "X-001",
    "service_code": "unknown_service",
})
assert resp.status_code == 400
print(f"  [OK] POST /order/unknown_service/submit -> 400 (expected)")

# 10. Non-existent order
resp = client.get("/order/NONEXISTENT/status")
assert resp.status_code == 404
print(f"  [OK] GET /order/NONEXISTENT/status -> 404 (expected)")

print(f"\n✅ All {10} API tests passed!")
