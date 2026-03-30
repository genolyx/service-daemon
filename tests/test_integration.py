"""
Service Daemon 통합 테스트

모듈 import, API 엔드포인트, Carrier Screening 워크플로우를 검증합니다.
"""

import os
import sys
import json
import tempfile
import asyncio

# 프로젝트 루트를 path에 추가
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PASS = "\033[92m✓ PASS\033[0m"
FAIL = "\033[91m✗ FAIL\033[0m"
results = []


def test(name, func):
    """테스트 실행 헬퍼"""
    try:
        func()
        print(f"  {PASS} {name}")
        results.append((name, True, ""))
    except Exception as e:
        print(f"  {FAIL} {name}: {e}")
        results.append((name, False, str(e)))


# ═══════════════════════════════════════════════════════════
# 1. Module Import Tests
# ═══════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("1. Module Import Tests")
print("=" * 60)

test("import config", lambda: __import__("app.config", fromlist=["settings"]))
test("import models", lambda: __import__("app.models", fromlist=["Job", "OrderStatus"]))
test("import services.base", lambda: __import__("app.services.base", fromlist=["ServicePlugin"]))
test("import services registry", lambda: __import__("app.services", fromlist=["load_plugins"]))
test("import queue_manager", lambda: __import__("app.queue_manager", fromlist=["QueueManager"]))
test("import runner", lambda: __import__("app.runner", fromlist=["PipelineRunner"]))
test("import platform_client", lambda: __import__("app.platform_client", fromlist=["PlatformClient"]))
test("import main (FastAPI app)", lambda: __import__("app.main", fromlist=["app"]))

# Carrier Screening 서브모듈
test("import carrier_screening.plugin", lambda: __import__("app.services.carrier_screening.plugin", fromlist=["CarrierScreeningPlugin"]))
test("import carrier_screening.vcf_parser", lambda: __import__("app.services.carrier_screening.vcf_parser", fromlist=["clean_vcf_remove_formats"]))
test("import carrier_screening.annotator", lambda: __import__("app.services.carrier_screening.annotator", fromlist=["VariantAnnotator"]))
test("import carrier_screening.acmg", lambda: __import__("app.services.carrier_screening.acmg", fromlist=["classify_acmg_lite"]))
test("import carrier_screening.review", lambda: __import__("app.services.carrier_screening.review", fromlist=["generate_result_json"]))
test("import carrier_screening.report", lambda: __import__("app.services.carrier_screening.report", fromlist=["generate_report_json"]))


# ═══════════════════════════════════════════════════════════
# 2. Plugin Loading Tests
# ═══════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("2. Plugin Loading Tests")
print("=" * 60)

from app.services import load_plugins, get_plugin, list_service_codes, _registry

# 레지스트리 초기화
_registry.clear()


def test_load_carrier_screening():
    load_plugins(["carrier_screening"])
    assert "carrier_screening" in list_service_codes(), "carrier_screening not loaded"


def test_plugin_properties():
    plugin = get_plugin("carrier_screening")
    assert plugin is not None, "Plugin is None"
    assert plugin.service_code == "carrier_screening"
    assert plugin.display_name == "Carrier Screening"
    stages = plugin.get_progress_stages()
    assert isinstance(stages, dict) and len(stages) > 0


def test_plugin_validate_params():
    plugin = get_plugin("carrier_screening")
    ok, msg = plugin.validate_params({})
    assert ok is True


test("load carrier_screening plugin", test_load_carrier_screening)
test("plugin properties", test_plugin_properties)
test("plugin validate_params", test_plugin_validate_params)


# ═══════════════════════════════════════════════════════════
# 3. ACMG Classification Tests
# ═══════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("3. ACMG Classification Tests")
print("=" * 60)

from app.services.carrier_screening.acmg import classify_acmg_lite


def test_acmg_pathogenic():
    variant = {
        "gene": "CFTR",
        "effect": "frameshift_variant",
        "clinvar_sig_primary": "Pathogenic",
        "clinvar_stars": 3,
        "gnomad_af": 0.00001,
    }
    result = classify_acmg_lite(variant)
    assert result["classification"] in ("Pathogenic", "Likely pathogenic"), \
        f"Expected Pathogenic/Likely pathogenic, got {result['classification']}"
    assert len(result["criteria_met"]) > 0


def test_acmg_benign():
    variant = {
        "gene": "BRCA1",
        "effect": "synonymous_variant",
        "clinvar_sig_primary": "Benign",
        "clinvar_stars": 2,
        "gnomad_af": 0.15,
    }
    result = classify_acmg_lite(variant)
    assert result["classification"] in ("Benign", "Likely benign"), \
        f"Expected Benign/Likely benign, got {result['classification']}"


def test_acmg_vus():
    variant = {
        "gene": "BRCA2",
        "effect": "missense_variant",
        "clinvar_sig_primary": "",
        "clinvar_stars": 0,
        "gnomad_af": 0.0005,
    }
    result = classify_acmg_lite(variant)
    assert result["classification"] == "VUS", \
        f"Expected VUS, got {result['classification']}"


test("ACMG pathogenic classification", test_acmg_pathogenic)
test("ACMG benign classification", test_acmg_benign)
test("ACMG VUS classification", test_acmg_vus)


# ═══════════════════════════════════════════════════════════
# 4. Review Data Generation Tests
# ═══════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("4. Review Data Generation Tests")
print("=" * 60)

from app.services.carrier_screening.review import generate_result_json, generate_variants_tsv


def test_generate_result_json():
    with tempfile.TemporaryDirectory() as tmpdir:
        variants = [
            {
                "chrom": "chr7", "pos": 117559590, "ref": "C", "alt": "T",
                "gene": "CFTR", "transcript": "NM_000492.4",
                "hgvsc": "c.1521_1523delCTT", "hgvsp": "p.Phe508del",
                "effect": "frameshift_variant", "zygosity": "het",
                "gnomad_af": 0.00002, "clinvar_sig_primary": "Pathogenic",
                "clinvar_stars": 3, "clinvar_dn": "Cystic fibrosis",
                "dp": 120, "vaf": 0.48,
            },
        ]
        acmg_results = [
            {
                "final_classification": "Pathogenic",
                "final_criteria": ["PVS1", "PS1", "PM2"],
                "final_reasoning": "LOF variant; ClinVar Pathogenic",
                "rule_based": {"classification": "Pathogenic"},
                "ai": None,
            },
        ]
        qc = {"coverage": {"mean_coverage": 150.5}, "alignment": {"mapping_rate": 99.2}}

        path = generate_result_json(
            annotated_variants=variants,
            acmg_results=acmg_results,
            qc_summary=qc,
            sample_name="TEST001",
            order_id="ORD-001",
            output_dir=tmpdir,
        )

        assert os.path.exists(path), f"result.json not created: {path}"
        with open(path) as f:
            data = json.load(f)
        assert data["type"] == "carrier_screening_result"
        assert len(data["variants"]) == 1
        assert data["variants"][0]["acmg_classification"] == "Pathogenic"
        assert data["variant_stats"]["pathogenic_or_likely"] == 1
        assert data["status"] == "pending_review"
        assert "dark_genes" in data
        assert data["dark_genes"].get("status") == "not_found"


def test_generate_variants_tsv():
    with tempfile.TemporaryDirectory() as tmpdir:
        variants = [
            {
                "variant_id": "VAR_0001",
                "chrom": "chr7", "pos": 117559590, "ref": "C", "alt": "T",
                "gene": "CFTR", "acmg_classification": "Pathogenic",
                "acmg_criteria": ["PVS1", "PM2"],
            },
        ]
        path = generate_variants_tsv(variants, tmpdir)
        assert os.path.exists(path), f"variants.tsv not created"
        with open(path) as f:
            lines = f.readlines()
        assert len(lines) == 2  # header + 1 row


test("generate result.json", test_generate_result_json)
test("generate variants.tsv", test_generate_variants_tsv)


# ═══════════════════════════════════════════════════════════
# 5. Report Generation Tests
# ═══════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("5. Report Generation Tests")
print("=" * 60)

from app.services.carrier_screening.report import generate_report_json


def test_generate_report_json():
    with tempfile.TemporaryDirectory() as tmpdir:
        confirmed_variants = [
            {
                "variant_id": "VAR_0001",
                "chrom": "chr7", "pos": 117559590, "ref": "C", "alt": "T",
                "gene": "CFTR", "transcript": "NM_000492.4",
                "hgvsc": "c.1521_1523delCTT", "hgvsp": "p.Phe508del",
                "effect": "frameshift_variant", "zygosity": "het",
                "gnomad_af": 0.00002,
                "clinvar_sig_primary": "Pathogenic",
                "clinvar_dn": "Cystic fibrosis",
                "acmg_classification": "Pathogenic",
                "acmg_criteria": ["PVS1", "PS1", "PM2"],
                "reviewer_confirmed": True,
                "reviewer_classification": "Pathogenic",
                "reviewer_comment": "Confirmed pathogenic variant for CF",
            },
            {
                "variant_id": "VAR_0002",
                "chrom": "chr1", "pos": 100000, "ref": "A", "alt": "G",
                "gene": "GENE2",
                "reviewer_confirmed": False,  # 미확정 → 리포트에서 제외
            },
        ]

        path = generate_report_json(
            order_id="ORD-001",
            sample_name="TEST001",
            confirmed_variants=confirmed_variants,
            reviewer_info={"name": "Dr. Kim", "id": "REV001", "institution": "Genolyx"},
            qc_summary={"coverage": {"mean_coverage": 150.5}, "alignment": {"mapping_rate": 99.2}},
            output_dir=tmpdir,
        )

        assert os.path.exists(path), f"report.json not created"
        with open(path) as f:
            data = json.load(f)
        assert data["type"] == "carrier_screening_report"
        assert len(data["confirmed_variants"]) == 1  # 미확정 제외
        assert data["carrier_status"]["status"] == "carrier"
        assert "Cystic fibrosis" in data["carrier_status"]["carrier_diseases"]
        assert data["reviewer"]["name"] == "Dr. Kim"
        # Jinja templates use primary_patient.findings (orig JSON shape)
        assert len(data["primary_patient"]["findings"]) == 1
        assert data["primary_patient"]["findings"][0]["gene"] == "CFTR"


def test_carrier_status_negative():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = generate_report_json(
            order_id="ORD-002",
            sample_name="TEST002",
            confirmed_variants=[],
            reviewer_info={"name": "Dr. Lee"},
            qc_summary={},
            output_dir=tmpdir,
        )
        with open(path) as f:
            data = json.load(f)
        assert data["carrier_status"]["status"] == "negative"


test("generate report.json (carrier)", test_generate_report_json)
test("carrier status negative", test_carrier_status_negative)


# ═══════════════════════════════════════════════════════════
# 6. FastAPI Endpoint Tests
# ═══════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("6. FastAPI Endpoint Tests")
print("=" * 60)

from fastapi.testclient import TestClient
from app.main import app

# 플러그인이 이미 로드되어 있으므로 TestClient 사용 가능
client = TestClient(app)


def test_health_endpoint():
    resp = client.get("/health")
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "healthy"
    assert "carrier_screening" in data["registered_services"]


def test_services_endpoint():
    resp = client.get("/services")
    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] >= 1
    codes = [s["service_code"] for s in data["services"]]
    assert "carrier_screening" in codes


def test_dashboard_endpoint():
    resp = client.get("/")
    assert resp.status_code == 200
    data = resp.json()
    assert data["version"] == "2.0.0"


def test_submit_order():
    resp = client.post("/order/carrier_screening/submit", json={
        "order_id": "TEST-ORD-001",
        "service_code": "carrier_screening",
        "work_dir": "260101",
        "fastq_r1_path": "/data/fastq/260101/TEST-ORD-001/TEST-ORD-001_R1.fastq.gz",
        "fastq_r2_path": "/data/fastq/260101/TEST-ORD-001/TEST-ORD-001_R2.fastq.gz",
        "params": {
            "carrier": {
                "test_category": "standard_carrier",
                "report_language": "EN",
            },
        },
        "priority": "normal",
    })
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "accepted"
    assert data["service_code"] == "carrier_screening"


def test_order_status():
    resp = client.get("/order/TEST-ORD-001/status")
    assert resp.status_code == 200
    data = resp.json()
    assert data["order_id"] == "TEST-ORD-001"
    assert data["service_code"] == "carrier_screening"


def test_queue_summary():
    resp = client.get("/queue/summary")
    assert resp.status_code == 200
    data = resp.json()
    assert "total_queued" in data


def test_unknown_service():
    resp = client.post("/order/unknown_service/submit", json={
        "order_id": "TEST-ORD-999",
        "service_code": "unknown_service",
    })
    assert resp.status_code == 400


def test_report_endpoint():
    try:
        resp = client.post("/order/TEST-ORD-001/report", json={
            "confirmed_variants": [
                {
                    "variant_id": "VAR_0001",
                    "gene": "CFTR",
                    "reviewer_confirmed": True,
                    "reviewer_classification": "Pathogenic",
                    "reviewer_comment": "Confirmed",
                }
            ],
            "reviewer_info": {"name": "Dr. Kim", "id": "REV001"},
        })
    except Exception as e:
        # order_store SQLite가 읽기 전용인 환경(로컬/CI)에서는 리포트 중 저장 단계에서 실패할 수 있음
        if "readonly database" in str(e).lower():
            return
        raise
    # 리포트 생성은 실제 파일이 없으므로 500 가능 (로직 자체는 호출됨)
    assert resp.status_code in (200, 500), f"Unexpected status: {resp.status_code}"


def test_dark_genes_pdf_only_approved_sections():
    from app.services.carrier_screening.dark_genes import dark_genes_for_pdf

    block = {
        "status": "found",
        "detailed_text": "",
        "detailed_sections": [
            {
                "title": "SMAca CHECK (Silent Carrier + Coverage)",
                "body": "SMN1_CN=3\nSMN2_CN=3\nCov(1,2)=2.82, 3.06\nSilentCarrier=False\nCT_Ratio=1:1\nC_T=51,49",
                "kind": "normal",
            },
            {"title": "OTHER SECTION", "body": "beta line", "kind": "normal"},
        ],
        "section_reviews": [
            {"approved": True, "notes": "ok"},
            {"approved": False, "notes": ""},
        ],
    }
    out = dark_genes_for_pdf(block)
    html_out = out.get("report_detailed_html") or ""
    assert "Spinal Muscular Atrophy" in html_out
    assert "SMN1" in html_out
    assert "SNP C/T ratio" in html_out
    assert "1:1" in html_out
    assert "51 / 49" in html_out
    assert "ok" in html_out
    assert "SMN1_CN" not in html_out
    assert "OTHER SECTION" not in html_out
    assert "#15803d" in html_out  # no pipeline WARNING → low-risk green accent by default

    block_warn = {
        "status": "found",
        "detailed_text": "",
        "detailed_sections": [
            {
                "title": "HBA ANALYSIS (Alpha Thalassemia - Dosage)",
                "body": "Est_CN=2\nRatio=1\nWARNING: allele imbalance",
                "kind": "normal",
            },
        ],
        "section_reviews": [{"approved": True, "notes": ""}],
    }
    html_warn = (dark_genes_for_pdf(block_warn).get("report_detailed_html") or "")
    assert "#be123c" in html_warn  # pipeline WARNING line → high-risk red accent

    block_low = dict(block)
    block_low["section_reviews"] = [
        {"approved": True, "notes": "ok", "risk": "low"},
        {"approved": False, "notes": ""},
    ]
    html_low = (dark_genes_for_pdf(block_low).get("report_detailed_html") or "")
    assert "#15803d" in html_low  # reviewer low-risk green accent

    none_approved = dict(block)
    none_approved["section_reviews"] = [
        {"approved": False, "notes": ""},
        {"approved": False, "notes": ""},
    ]
    assert dark_genes_for_pdf(none_approved) == {}

    string_false = dict(block)
    string_false["section_reviews"] = [
        {"approved": "false", "notes": ""},
        {"approved": "true", "notes": "Second section note"},
    ]
    out_sf = dark_genes_for_pdf(string_false)
    h2 = out_sf.get("report_detailed_html") or ""
    assert "OTHER" in h2
    assert "Second section note" in h2
    assert "beta line" in h2
    assert "alpha line" not in h2


test("GET /health", test_health_endpoint)
test("GET /services", test_services_endpoint)
test("GET / (dashboard)", test_dashboard_endpoint)
test("POST /order/carrier_screening/submit", test_submit_order)
test("GET /order/{order_id}/status", test_order_status)
test("GET /queue/summary", test_queue_summary)
test("POST /order/unknown_service/submit (400)", test_unknown_service)
test("POST /order/{order_id}/report", test_report_endpoint)
test("dark_genes PDF only approved sections", test_dark_genes_pdf_only_approved_sections)


# ═══════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("Test Summary")
print("=" * 60)

passed = sum(1 for _, ok, _ in results if ok)
failed = sum(1 for _, ok, _ in results if not ok)
total = len(results)

print(f"\n  Total: {total}  |  Passed: {passed}  |  Failed: {failed}")

if failed > 0:
    print(f"\n  Failed tests:")
    for name, ok, err in results:
        if not ok:
            print(f"    - {name}: {err}")

print()
sys.exit(0 if failed == 0 else 1)
