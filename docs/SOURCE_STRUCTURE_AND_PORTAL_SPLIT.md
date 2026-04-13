# Service-Daemon 소스 구조 및 Portal 분리 가이드

이 문서는 **service-daemon** 레포지토리의 디렉터리·책임을 정리하고, **기존 Portal에 기능을 얹을 때** 무엇을 데몬에 두고 무엇을 UI 쪽으로 가져가야 하는지 구분하기 위한 참고 자료입니다.

---

## 1. 한눈에 보는 구조

이 프로젝트는 **FastAPI 백엔드**와 **단일 HTML Portal(`portal/index.html`)**이 **같은 프로세스**에서 서빙되는 형태입니다.

| 층 | 위치 | 역할 |
|----|------|------|
| **Portal UI** | `portal/index.html` (단일 파일, JS 인라인) | 주문·큐·Review·문헌 등 브라우저 UI. API를 호출만 함. |
| **HTTP API** | `app/main.py` (+ 일부 라우트 보조 모듈) | 주문/큐/결과/파일/리포트/리뷰 저장 등 REST 엔드포인트. |
| **오케스트레이션** | `app/queue_manager.py`, `app/runner.py` | 작업 큐, 동시 실행 제한, 파이프라인 실행 루프. |
| **서비스별 로직** | `app/services/<service_code>/` | 플러그인: 입력 준비, 파이프라인 명령, 완료 판정, `process_results`. |
| **공통 도메인 로직** | `app/services/carrier_screening/*.py` 등 | VCF 처리, `result.json` 생성, 리포트, Dark genes, PGx 등. |
| **설정·모델** | `app/config.py`, `app/models.py` | 환경 변수, Pydantic 요청/응답 모델. |
| **데이터 자산** | `data/report_templates/`, `data/db/` 등 | PDF HTML 템플릿, TSV 등 (런타임 읽기). |
| **외부 파이프라인** | 이 레포 **밖** (예: carrier-screening Nextflow) | 실제 시퀀싱 파이프라인; 데몬은 명령 실행·산출물 경로만 다룸. |

**중요:** Review 탭에 보이는 변이·QC·Coverage 등 **데이터의 “진실”은 대부분 서버가 만든 `result.json`**이며, Portal은 이를 **GET `/order/{id}/result`** 로 받아 표시합니다. UI만 옮기고 API를 붙이지 않으면 동일 기능을 재현할 수 없습니다.

---

## 2. 디렉터리 트리 (요약)

```
service-daemon/
├── app/
│   ├── main.py                 # FastAPI 앱, 거의 모든 라우트
│   ├── config.py               # Settings (pydantic-settings, .env)
│   ├── models.py               # 요청/응답·Job 등 Pydantic 모델
│   ├── queue_manager.py        # 큐·워커·Job 상태
│   ├── runner.py               # 파이프라인 서브프로세스 실행
│   ├── order_store.py          # 주문 영속화·리포트 JSON 동기화 등
│   ├── order_cleanup.py        # 삭제/정리
│   ├── platform_client.py      # Genolyx Platform API (선택)
│   ├── auth_client.py          # 인증 (선택)
│   ├── annotation_resources.py # 주석 리소스 점검 (/health?annotation=true)
│   ├── logging_config.py
│   ├── datetime_kst.py
│   ├── telegram_notify.py      # (선택) 텔레그램 알림
│   ├── services/
│   │   ├── __init__.py         # 플러그인 로드·레지스트리
│   │   ├── base.py             # ServicePlugin 추상 클래스
│   │   ├── carrier_screening/  # Carrier / WES / Health screening 공통 스택
│   │   ├── whole_exome/        # 플러그인 래퍼
│   │   ├── health_screening/
│   │   ├── sgnipt.py
│   │   ├── nipt.py
│   │   ├── wes_panels.py       # 패널 카탈로그·해석 유전자 집합
│   │   ├── gene_panel_coverage.py
│   │   └── _resolve_ai.py
│   └── static/                 # (있다면) 기타 정적 파일
├── portal/
│   └── index.html              # 전체 Portal UI (대용량 단일 파일)
├── data/
│   ├── report_templates/       # WeasyPrint용 HTML (carrier, exome, pgx …)
│   ├── db/                     # 예: PGx 커스텀 TSV 등
│   └── …
├── docs/                       # 문서
├── scripts/                    # 운영·동기화 스크립트 (선택)
├── tests/                      # pytest
├── docker-compose.yml
├── docker-compose-dev.yml
├── Dockerfile
├── requirements.txt
├── run-integration.sh
├── README.md
└── ARCHITECTURE.md             # 플러그인·큐 개념 (다소 초기 버전)
```

---

## 3. 실행 시 데이터 흐름 (Review 기능 기준)

1. **파이프라인**이 샘플 디렉터리에 VCF·BAM 등을 남김 (외부 Nextflow 등).
2. 플러그인 **`process_results`** (`carrier_screening/plugin.py` 등)가 VCF를 읽고 주석·ACMG·QC를 합쳐 **`output/.../result.json`** 을 생성 (`carrier_screening/review.py` 등).
3. Portal **Load** → **`GET /order/{order_id}/result`** → 디스크의 `result.json`을 JSON으로 반환 (필요 시 후처리·dark_genes 보강).
4. **Coverage 탭** → BAM은 **`GET /order/{id}/file/...`** 로 범위 요청, 유전자 커버리지는 **`/gene-coverage/{gene}`**, **`/coverage-context`** 등.
5. **리뷰어 입력 저장** → **`PUT/PATCH .../gene-knowledge`**, **`variant-knowledge`**, **`dark-genes-review`**, **`pgx-review`** 등이 서버 SQLite/파일에 반영.

즉 **“가져갈 UI”는 `portal/index.html` 안의 JS·HTML**이지만, **동일 동작을 원하면 동일 API와 동일 `result.json` 스키마**를 소비해야 합니다.

---

## 4. `app/main.py` 라우트 그룹 (참고용)

파일이 크므로, **기능별로 묶인 엔드포인트**만 기억하면 통합 설계에 유리합니다.

| 그룹 | 예시 경로 | 용도 |
|------|-----------|------|
| 헬스 | `GET /health` | 상태·등록 서비스·FASTQ 루트 |
| 주문 CRUD·실행 | `POST .../save`, `PATCH /order/{id}`, `POST .../submit`, `POST .../start`, `POST .../stop`, `DELETE ...` | 주문 생명주기 |
| 결과·리뷰 | `GET /order/{id}/result`, `POST .../reprocess-results` | `result.json` 조회·재생성 |
| 파일 | `GET /order/{id}/files`, `GET .../file/{path}` | BAM·VCF·PDF 등 |
| Review 보조 | `GET .../gene-coverage/...`, `GET .../coverage-context` | Coverage 탭 |
| 리뷰 저장 | `PUT .../gene-knowledge`, `PUT .../variant-knowledge`, `POST/PATCH .../dark-genes-review`, `.../pgx-review` | Portal에서 편집한 메타데이터 |
| 리포트 | `POST .../report`, `.../report/preview`, `.../from-html` | PDF/HTML 생성 |
| Portal 보조 API | `GET /api/portal/wes-panels`, `GET /api/portal/resources`, FASTQ browse 등 | 주문 폼·리소스 |
| 문헌 | `GET /api/literature/*` | Literature 뷰 |
| 큐 | `GET /queue/summary`, `GET /queue/status` | 대시보드 |

전체 목록은 `app/main.py`에서 `@app.get` / `@app.post` 등으로 검색하는 것이 가장 정확합니다.

---

## 5. 서비스 플러그인 (`app/services/`)

- **`base.ServicePlugin`**: 모든 서비스가 구현해야 하는 메서드 (`prepare_inputs`, `get_pipeline_command`, `check_completion`, `process_results`, `get_output_files`, `generate_report` 등).
- **`load_plugins(enabled_services)`**: `ENABLED_SERVICES` 환경 변수에 따라 `app.services.<code>` 모듈을 로드하고 `create_plugin()`으로 인스턴스 등록.
- **`carrier_screening/`**: Carrier·Whole exome·Health screening 계열의 **핵심 비즈니스 로직** (VCF, 주석, `result.json`, 리포트, dark genes, PGx, 문헌 연동 등).
- **`sgnipt.py`**, **`nipt.py`**: 서비스별 `result.json` 규칙이 다름.

새 서비스를 **데몬 안에만** 추가할 때는 이 플러그인 패턴을 따릅니다. **외부 Portal**에서는 이 Python 코드를 “가져가기”보다 **HTTP로 호출**하는 것이 맞습니다.

---

## 6. 데몬에만 두는 것 vs Portal(또는 프론트)로 가져가는 것

### 6.1 반드시 데몬(백엔드)에 남겨야 하는 것

| 영역 | 이유 |
|------|------|
| `queue_manager`, `runner`, 플러그인 `process_results` | 로컬 경로·서브프로세스·동시성·장시간 작업은 브라우저에서 불가능. |
| `result.json` 생성·재처리 (`reprocess-results`) | VCF/pysam/주석 DB 접근·대용량 I/O. |
| `GET /order/.../file/...` (Range 요청 등) | BAM 스트리밍·보안·경로 해제. |
| 리포트 생성 (`report.py` + `data/report_templates/`) | WeasyPrint, Jinja, 파일 쓰기. |
| `gene_knowledge` / `variant_knowledge` SQLite·저장 로직 | 서버 측 영속화·권한. |
| Platform 업로드·인증 (`platform_client`, `api_key`) | 시크릿·서버 간 통신. |

### 6.2 “가져가도 되는” 것 (UI/클라이언트 측)

| 영역 | 설명 |
|------|------|
| **`portal/index.html`의 HTML/CSS/JS** | 주문 폼, Orders 테이블, Review 탭, Literature 등 **표현층**. 다른 SPA/프레임워크로 옮길 때 **참고·포팅** 대상. |
| **API 호출 규약** | `API_BASE`, `ORDER_API_PREFIX`, `apiGet`/`apiPost` 패턴 — 새 Portal에서 **동일 URL**로 호출하면 동작을 재현하기 쉬움. |
| **표시 전용 로직** | 테이블 정렬, IGV.js 초기화, 탭 전환 등 **순수 UI** (단, BAM URL은 데몬이 준 `file` URL을 써야 함). |

### 6.3 경계에 있는 것 (팀 합의 필요)

| 영역 | 비고 |
|------|------|
| **`result.json` 스키마** | UI는 스키마에 강하게 결합되어 있음. 백엜을 그대로 쓰면 **`carrier_screening/review.py` 문서화 주석**과 실제 JSON을 **계약(contract)** 으로 보는 것이 좋음. |
| **Literature (`/api/literature/*`)** | 데몬에 구현되어 있으면 프록시만 해도 되고, 별도 마이크로서비스로 쪼갤 수도 있음. |
| **WES 패널 API** | `wes_panels.py` + 라우트 — UI는 카탈로그 JSON만 알면 됨. |

---

## 7. 기존 Portal에 기능을 “얹는” 통합 패턴

1. **API 전용 통합 (권장)**  
   - 데몬은 지금처럼 **별도 origin** (`https://daemon.example.com`)으로 두고, 기존 Portal 프론트에서 **CORS** 또는 **같은 도메인 리버스 프록시**(`/api/order/...` → 데몬)로 연결.  
   - 이 레포의 `ORDER_API_PREFIX`·`apiPath()` 설계가 프록시 경로를 고려함 (`portal/index.html` 상단 주석 참고).

2. **iframe / 링크**  
   - 빠르게 통합할 때: 데몬이 서빙하는 `/portal/` 을 **새 창 또는 iframe**으로 임베드. 스타일·SSO는 별도 과제.

3. **UI만 복사**  
   - `portal/index.html`에서 필요한 **섹션·함수**만 추출해 기존 앱으로 이식.  
   - 반드시 **같은 API 베이스 URL**과 **인증 헤더**(`Authorization` / `X-API-Key` — `config.py`의 `api_key`)를 맞출 것.

---

## 8. 파일별 “참고 우선순위” (포팅 체크리스트)

| 목표 | 우선 볼 파일 |
|------|----------------|
| 주문·큐·상태 | `app/main.py` (주문 관련), `app/queue_manager.py`, `app/models.py` |
| Review 데이터 생성 | `app/services/carrier_screening/review.py`, `plugin.py` (`process_results`) |
| Review UI 동작 | `portal/index.html` (`loadReviewData`, `showReviewTab`, `renderReviewQC`, Coverage/IGV 부분) |
| Coverage API | `app/services/gene_panel_coverage.py`, `main.py`의 `gene-coverage`, `coverage-context` |
| Dark genes | `app/services/carrier_screening/dark_genes.py` + Portal 내 dark genes 렌더 |
| PGx | `app/services/carrier_screening/pgx_report.py` + `pgx-review` 라우트 |
| 리포트 PDF | `app/services/carrier_screening/report.py`, `data/report_templates/` |
| sgNIPT | `app/services/sgnipt.py` |
| 설정 예시 | `.env.example`, `README.md` |

---

## 9. 관련 문서

- `ARCHITECTURE.md` — 플러그인·큐 개념 (초기 설계; 디렉터리 이름은 현재 트리와 다를 수 있음).
- `README.md` — 연동 테스트·포트·Portal URL.
- `docs/integrations/carrier_dark_genes.md` — Dark genes·Portal 연동 세부.

---

## 10. 요약

| 질문 | 답 |
|------|-----|
| Portal 기능만 다른 저장소로 옮기려면? | **`portal/index.html`** 을 참고해 UI/JS를 옮기고, **데몬 API URL**을 그대로 쓴다. |
| “로직”은 어디 있나? | **대부분 `app/services/` + `app/main.py`**. `portal/`에는 비즈니스 로직이 거의 없다. |
| 데몬 없이 Portal만으로 될까? | **동일 기능 불가**. `result.json` 생성·파일 서빙·리포트는 서버 필요. |
| 가장 중요한 계약은? | **`GET /order/{id}/result` 의 JSON 스키마**와 주문·파일 관련 API 목록. |

이 문서는 레포 구조가 바뀔 수 있으므로, **라우트와 플러그인 목록은 `app/main.py`와 `app/services/`를 최종 기준**으로 확인하세요.
