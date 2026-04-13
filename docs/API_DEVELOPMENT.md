# Service Daemon — Portal 연동·개발용 API 안내

데몬을 띄운 뒤 **Portal UI**와 **API 문서**를 같이 쓰면서 개발할 때 참고합니다.

## 정본은 HTML (Swagger / ReDoc)

**경로·파라미터·요청/응답 스키마가 정확히 맞아야 할 때는 이 레포의 Markdown이 아니라, 데몬이 띄우는 HTML 문서를 쓰는 것이 맞습니다.** FastAPI가 코드에서 **OpenAPI**를 생성하고, 그걸 렌더링한 것이 아래입니다.

| 형식 | URL | 설명 |
|------|-----|------|
| **Swagger UI** | `http://localhost:<APP_PORT>/docs` | 메서드별 **Try it out**, 스키마 전개 |
| **ReDoc** | `http://localhost:<APP_PORT>/redoc` | 읽기 좋은 한 페이지 레퍼런스 |
| **OpenAPI JSON** | `http://localhost:<APP_PORT>/openapi.json` | 기계 처리·다른 도구 import |

브라우저로 **`/docs`** 만 열어도 됩니다. 이 문서(`.md`)는 **워크플로 설명·Portal이 쓰는 API 요약**용이며, **스펙의 단일 진실 공급원(SSOT)은 OpenAPI**입니다.

`curl`로 루트만 조회할 때도 링크가 나옵니다: `curl -s -H 'Accept: application/json' http://localhost:<APP_PORT>/` → `api_docs` 필드.

## 빠른 링크 (포트는 환경에 맞게 바꿈)

| 용도 | URL |
|------|-----|
| **Portal** | `http://localhost:<APP_PORT>/portal/` |
| **Swagger UI** | `http://localhost:<APP_PORT>/docs` |
| **ReDoc** | `http://localhost:<APP_PORT>/redoc` |
| **OpenAPI JSON** | `http://localhost:<APP_PORT>/openapi.json` |
| **헬스** | `http://localhost:<APP_PORT>/health` |

로컬 실행 예: `./run-integration.sh` → 터미널에 찍힌 포트 사용.

## Portal이 실제로 호출하는 API (전체 목록이 아님)

`docs/API_DEVELOPMENT.md` 아래쪽 **엔드포인트 요약**은 데몬이 노출하는 **전체에 가깝고**, 내장 **`portal/index.html`이 쓰는 것은 그중 일부**입니다. 다른 클라이언트(플랫폼, 스크립트, 자동화)는 추가 엔드포인트를 쓸 수 있습니다.

아래는 **`portal/index.html` 기준** (같은 오리진에서 `fetch` / `apiGet` 등으로 호출하는 경로)입니다.

### 대시·헬스

- `GET /health`
- `GET /services`
- `GET /queue/summary`
- `GET /queue/dashboard-bucket` (쿼리: `bucket`, `sort`, `order`, 선택 `service_code`)

### 주문 (저장·실행·정리)

- `POST /order/{service_code}/save`
- `PATCH /order/{order_id}` (주문 본문 수정)
- `POST /order/{order_id}/start` (본 UI의 **Submit / Force Run** 등 — 바디에 `fresh` 등)
- `POST /order/{order_id}/reprocess-results`
- `POST /order/{order_id}/stop`
- `POST /order/{order_id}/delete-run`
- `POST /order/{order_id}/purge-db` (쿼리 있을 수 있음)
- `GET /orders`
- `GET /order/{order_id}`

### Review·파일·리포트

- `GET /order/{order_id}/result` (`debug=1` 포함 가능)
- `GET /order/{order_id}/files`
- `GET /order/{order_id}/file/{filename:path}` (직접 `fetch`로 URL 조합 — BAM 등)
- `GET /order/{order_id}/pipeline-log` (텍스트)
- `GET /order/{order_id}/gene-coverage/{gene_symbol}`
- `GET /order/{order_id}/coverage-context`
- `GET /order/{order_id}/gene-knowledge` (쿼리: `gene`, `force` 등)
- `PUT /order/{order_id}/gene-knowledge`
- `PUT /order/{order_id}/variant-knowledge`
- `POST /order/{order_id}/dark-genes-review`
- `POST /order/{order_id}/pgx-review`
- `POST /order/{order_id}/report`
- `POST /order/{order_id}/report/preview`
- `POST /order/{order_id}/report/from-html`

### Portal 보조·패널·리소스

- `GET /api/portal/wes-panels`
- `POST /api/portal/wes-panels/custom`
- `DELETE /api/portal/wes-panels/custom/{panel_id}`
- `GET /api/portal/resources`

### FASTQ / BAM CSV 브라우즈

- `GET /api/fastq/browse` (실패 시 `GET /api/portal/browse/fastq` 재시도)
- `GET /api/portal/bam-csv/browse` (또는 `abs_path` 등)
- `GET /api/portal/bam-csv/sample-ids` (`fetch`로 직접 호출)

### Literature 탭

- `GET /api/literature/search`
- `GET /api/literature/articles`
- `GET /api/literature/articles/stats`
- `GET /api/literature/articles/{pmid}`

### 이 Portal에서 안 쓰는 것 (예시)

다음은 OpenAPI에는 있으나 **현재 `portal/index.html`에서는 호출하지 않는** 대표 예입니다. 통합 시 필요하면 별도로 붙입니다.

- `POST /order/{service_code}/submit` — UI는 **저장(`save`) 후 `start`** 조합을 사용
- `GET /order/{order_id}/status`, `POST /order/{order_id}/cancel` — 상태는 **`GET /orders`** 등으로 갱신
- `PATCH /order/{order_id}/fastq` — FASTQ는 주로 **저장/PATCH 주문**에 포함
- `GET /queue/status`
- `GET /api/bam-csv/browse` (별칭) — UI는 **`/api/portal/bam-csv/browse`** 우선
- Literature 캐시 삭제 `DELETE /api/literature/...` 등
- `POST /test/inject-mock-job` — 개발·테스트용

## 인증 (`API_KEY`)

`.env`에 **`API_KEY`** 가 비어 있으면 모든 API를 인증 없이 호출할 수 있습니다 (개발에 편함).

값이 설정되어 있으면 대부분의 엔드포인트에 다음 중 하나가 필요합니다.

```http
Authorization: Bearer <API_KEY>
```

또는

```http
X-API-Key: <API_KEY>
```

**키 없이 접근 가능한 예:** `/health`, `/portal/*`, **`/docs`**, **`/openapi.json`**, **`/redoc`**, 일부 `/api/portal/*`, `/static/*` 등. (전체는 `app/main.py`의 `api_access_key_guard` 참고.)

Swagger에서 **Authorize** 버튼에 `Bearer <API_KEY>` 또는 API Key 방식으로 넣고 호출할 수 있습니다.

## 리버스 프록시와 `/api` 접두사

앞단에서 **`/api/order/...`** 만 데몬으로 넘기는 경우, 서버는 동일 동작을 **`/order/...`** 로 다시 매핑합니다. OpenAPI에 적힌 경로는 기본적으로 **`/order/...`** 기준이므로, 프록시 뒤에서는 **`/api` 를 앞에 붙인 URL**로 호출하면 됩니다.

예: `GET /api/order/{order_id}/result` → 내부적으로 `GET /order/{order_id}/result`.

예외: `/api/literature/*`, `/api/portal/*`, `/api/fastq/*` 등은 원래부터 `/api` 가 경로의 일부입니다.

## 엔드포인트 요약 (그룹별)

아래는 자주 쓰는 순서로 묶었습니다. **상세 스키마·쿼리 파라미터는 `/docs` 가 정본**입니다.

### 헬스·메타

| Method | Path | 설명 |
|--------|------|------|
| GET | `/health` | 상태, 등록 서비스, FASTQ 루트 등 |
| GET | `/services` | 활성 플러그인 `service_code` 목록 |
| GET | `/` | 루트 (링크 안내 등) |

### 주문·실행

| Method | Path | 설명 |
|--------|------|------|
| POST | `/order/{service_code}/save` | 초안 저장 |
| PATCH | `/order/{order_id}` | 주문 수정 |
| POST | `/order/{service_code}/submit` | 제출 |
| POST | `/order/{order_id}/start` | 큐에서 실행 시작 |
| POST | `/order/{order_id}/stop` | 중지 |
| POST | `/order/{order_id}/cancel` | 취소 |
| POST | `/order/{order_id}/delete-run` | 런 삭제 |
| POST | `/order/{order_id}/purge-db` | DB 정리 |
| GET | `/order/{order_id}` | 주문 단건 |
| GET | `/orders` | 목록 |
| GET | `/order/{order_id}/status` | 상태 |
| PATCH | `/order/{order_id}/fastq` | FASTQ 경로 수정 |
| POST | `/order/{order_id}/reprocess-results` | 파이프라인 없이 결과만 재생성 |

### Review·결과

| Method | Path | 설명 |
|--------|------|------|
| GET | `/order/{order_id}/result` | `result.json` (Review Load) |
| GET | `/order/{order_id}/gene-coverage/{gene_symbol}` | 유전자 커버리지 |
| GET | `/order/{order_id}/coverage-context` | IGV·BAM 컨텍스트 |
| PUT | `/order/{order_id}/gene-knowledge` | Gene 지식 저장 |
| GET | `/order/{order_id}/gene-knowledge` | Gene 지식 조회 |
| PUT | `/order/{order_id}/variant-knowledge` | Variant 지식 저장 |
| PATCH/POST | `/order/{order_id}/dark-genes-review` | Dark genes 리뷰 |
| PATCH/POST | `/order/{order_id}/pgx-review` | PGx 리뷰 |

### 파일

| Method | Path | 설명 |
|--------|------|------|
| GET | `/order/{order_id}/files` | 파일 목록 |
| GET | `/order/{order_id}/file/{filename:path}` | 단일 파일 (BAM 등 Range 지원) |
| GET | `/order/{order_id}/pipeline-log` | 파이프라인 로그 텍스트 |

### 리포트

| Method | Path | 설명 |
|--------|------|------|
| POST | `/order/{order_id}/report` | 리포트 생성 |
| POST | `/order/{order_id}/report/preview` | HTML 프리뷰 |
| POST | `/order/{order_id}/report/from-html` | HTML 기반 생성 |

### 큐

| Method | Path | 설명 |
|--------|------|------|
| GET | `/queue/summary` | 큐 요약 |
| GET | `/queue/status` | 상태 |
| GET | `/queue/dashboard-bucket` | 대시보드 버킷 |

### Portal 보조

| Method | Path | 설명 |
|--------|------|------|
| GET | `/api/portal/wes-panels` | WES 패널 카탈로그 |
| POST | `/api/portal/wes-panels/custom` | 커스텀 패널 저장 |
| DELETE | `/api/portal/wes-panels/custom/{panel_id}` | 삭제 |
| GET | `/api/wes-panels` | 패널 (별칭) |
| GET | `/api/portal/resources` | 리소스 목록 |
| GET | `/api/portal/browse/fastq` 등 | FASTQ/BAM CSV 브라우즈 |

### 문헌 (Literature 탭)

| Method | Path |
|--------|------|
| GET | `/api/literature/search` |
| GET | `/api/literature/articles` |
| GET | `/api/literature/articles/stats` |
| GET | `/api/literature/articles/{pmid}` |
| DELETE | `/api/literature/articles/{pmid}` |
| DELETE | `/api/literature/cache` |
| GET | `/api/literature/searches` |

### FASTQ·BAM CSV (직접 경로)

| Method | Path |
|--------|------|
| GET | `/api/fastq/browse` |
| GET | `/api/bam-csv/browse` |
| GET | `/api/portal/bam-csv/browse` |
| GET | `/api/portal/bam-csv/sample-ids` |

### 테스트 전용

| Method | Path | 설명 |
|--------|------|------|
| POST | `/test/inject-mock-job` | 모의 잡 주입 (개발·테스트) |

## curl 예시

```bash
BASE=http://localhost:8003
# API_KEY가 있을 때
# curl -H "Authorization: Bearer $API_KEY" "$BASE/orders"

curl -sS "$BASE/health" | jq .
curl -sS "$BASE/services" | jq .
curl -sS "$BASE/orders" | jq .
```

## 소스 구조와의 관계

- **OpenAPI 스키마**는 코드의 타입 힌트·`response_model`에서 생성됩니다.
- Portal이 호출하는 경로와의 대응은 **`docs/SOURCE_STRUCTURE_AND_PORTAL_SPLIT.md`** 를 참고하세요.
