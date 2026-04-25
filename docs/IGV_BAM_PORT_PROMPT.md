# 다른 프로젝트용: BAM + IGV 기능 구현 프롬프트 (service-daemon 기준)

이 문서는 **service-daemon** 저장소에 이미 구현된 BAM·IGV 흐름을 바탕으로, 다른 프로젝트에 이식할 때 사용할 수 있는 **상세 개발 프롬프트**입니다.

이 저장소에는 **BAM을 통한 IGV**가 크게 두 갈래로 나뉩니다.

1. **브라우저에서 IGV.js로 BAM/BAI를 직접 로드**하는 Coverage 탭
2. **파이프라인이 만든 igv-reports HTML**을 가져와 iframe/모달로 보여주는 방식

---

## 배경 및 목표

- 주문(또는 케이스) 단위로 **인덱스된 BAM**(`*.bam` + 같은 경로의 `*.bai`)을 찾아, 웹 UI에서 **IGV.js**로 정렬(alignment) 트랙을 표시한다.
- 유전자/좌표로 점프할 수 있게 하고, **전체 염색체 줌**에서는 읽기가 비어 보일 수 있음을 UX로 안내한다 (IGV.js는 좁은 구간에서 pileup이 의미 있음).
- 선택적으로 **정적 igv-reports HTML** 리포트를 서버에서 fetch한 뒤 blob URL로 iframe에 넣는 흐름도 지원한다.

---

## 1. 백엔드: BAM 발견 및 `coverage-context` API

**목적**: 프론트가 “어떤 BAM을 어떤 URL로 로드할지”를 알 수 있게 JSON을 준다.

**엔드포인트 예시**: `GET /order/{order_id}/coverage-context`

### 응답 필드 (권장)

- `order_id`
- `genome_id`: 예 `"hg38"` (IGV.js `createBrowser`의 `genome`과 일치)
- `interpretation_genes`: 패널 해석에 쓰는 유전자 심볼 목록 (검색/필터용; 서비스별로 비울 수 있음)
- `bam_tracks`: 배열. 각 원소는 최소:
  - `rel_path`: 주문 산출물 루트 기준 **상대 경로** (슬래시 `/`, `..` 없음)
  - `label`: 표시용 파일명
  - `has_index`: `true` iff `*.bai`가 같은 루트 아래에 존재
  - `index_rel_path`: BAI의 상대 경로 (없으면 `null` 또는 생략)
  - (선택) `ancillary`: `true`면 “보조/특수 목적 BAM” (paraphase, SMN, FMR1 등)

### BAM 탐색 우선순위 (이 프로젝트 로직 요약)

1. **명시 경로**: 주문 파라미터 `params.igv_bam` 또는 `params.carrier.igv_bam`에 **절대 경로**가 있으면, 그 파일이 존재하고 주문 산출물 루트 아래로 매핑 가능하면 최우선 트랙으로 사용.
2. **mosdepth 형제**: `mosdepth` 산출물 `*.per-base.bed.gz` 같은 파일의 **stem**과 같은 디렉터리(상위 몇 단계까지)에서 stem이 맞는 `*.bam`을 찾고, 파일명 휴리스틱으로 “메인 exome BAM”을 고른다.
3. **재귀 glob**: 주문 관련 루트들 아래 `**/*.bam`을 스캔해 후보를 모은 뒤, 파일명·수정시간으로 정렬.

### 보조(ancillary) BAM 제외 규칙

- 파일명에 특정 토큰이 포함되면 (예: `paraphase`, `smn`, `fmr1`, `fmr1`, `eh_realigned`, `strc`, `gba`, `hba`, `no_dup` 등) **ancillary**로 표시.
- **자동 IGV용 목록**에서는 `ancillary`가 아닌 트랙만 노출하거나, ancillary만 있을 때 사용자에게 `igv_bam_message` 같은 문자열로 “메인 BAM을 `params.igv_bam` 또는 분석 출력 경로에 두라”고 안내.

### 인덱스 요구

IGV.js는 BAI가 필요하므로 `has_index === false`인 트랙은 자동 로드 후보에서 제외.

### prior reuse / 재분석

이전 주문 산출물을 재사용하는 경우, **이전 러닝의 analysis/output 트리**를 artifact 루트 목록 앞에 포함해 BAM·mosdepth·파일 경로가 동일한 정렬을 보게 한다.

---

## 2. 백엔드: 유전자 커버리지 → IGV locus

**목적**: 유전자 심볼로 조회하면 **BED 구간**을 돌려주고, 프론트가 `chrom:start-end` 형태로 IGV에 넘긴다.

**엔드포인트 예시**: `GET /order/{order_id}/gene-coverage/{gene_symbol}`

- 응답에 `bed_regions` (각각 `chrom`, `start`, `end`)가 있으면, 프론트는 첫 구간들을 합쳐  
  `chrom:minStart+1-maxEnd` (1-based closed interval 스타일) 또는 `gene_symbol` 문자열로 `search`한다.

---

## 3. 백엔드: BAM/BAI HTTP 제공 (필수 동작)

**엔드포인트 예시**: `GET /order/{order_id}/file/{filename:path}`  
(`filename`에 슬래시 포함 — `qc/foo.bam` 형태)

### 보안

- `path`에 `..` 금지, 실제 경로가 허용된 **주문 루트들**의 하위인지 `realpath`로 검증.
- 동일 상대 경로가 여러 루트에 있으면 **mtime 최신** 파일을 선택 (리포트와 일치).

### IGV.js 호환

- **HEAD** ` /order/{order_id}/file/...` 에서 `Content-Length`, `Content-Type`, **`Accept-Ranges: bytes`** 를 반환 (클라이언트가 길이를 먼저 조회).
- **GET**은 Range 요청을 처리할 수 있어야 BAM 스트리밍이 안정적이다 (Starlette/FastAPI의 `FileResponse` 등이 Range를 지원하는지 확인하고, 미지원이면 Range-aware 응답으로 교체).

### 인증

API 키가 있다면, 프론트에서 BAM/BAI fetch 시 **`X-API-Key`**(또는 동일한 정책)를 헤더로 붙이고, IGV 트랙 정의에 `headers`로 전달한다.

---

## 4. 프론트엔드: IGV.js (Coverage 탭 패턴)

**라이브러리**: CDN 예 — `igv@2.15.5` (`igv.min.js`), `window.igv` 로드 후 `igv.createBrowser` 사용.

**컨텍스트 로드**: Coverage 탭 진입 또는 주문 로드 시 `GET .../coverage-context`로 `reviewCoverageContext` 저장.

**트랙 선택**: `bam_tracks` 중 **`has_index === true`** 인 첫 트랙을 “자동 primary”로 사용.

### URL 생성

- `bamUrl` = `{API_BASE}/order/{orderId}/file/{encodeURIComponent 각 세그먼트를 / 로 join}`  
  (경로 전체를 한 번에 encode하지 말고 **세그먼트별** encode 권장)
- `indexURL` = 동일하게 `index_rel_path` 사용

### 트랙 객체

```json
{
  "type": "alignment",
  "format": "bam",
  "url": "<bamUrl>",
  "indexURL": "<indexUrl>",
  "name": "<label>",
  "height": "<동적>",
  "maxHeight": "<height + 여유>"
}
```

API 키가 있으면 `track.headers = { "X-API-Key": "..." }`.

### 브라우저 옵션

- `genome`: `coverage-context.genome_id` (예 `hg38`)
- `locus`: 초기 locus — 비어 있으면 **첫 변이 좌표 기준 좁은 창** (예 pos ± 45bp) 또는 유전자 심볼
- `tracks`: 위 alignment 트랙 하나

### 레이아웃

컨테이너를 `#coverageIgvHost`처럼 **밝은 배경**이고 충분한 높이(예 min-height)를 주고, `.igv-container`가 세로로 꽉 차게 CSS 조정.  
alignment 트랙 높이는 `ResizeObserver`로 호스트 높이에 맞춰 `setTrackHeight` 호출.

### 수명 관리

주문 전환 시 `browser.dispose()`, ResizeObserver `disconnect()`, 호스트 `innerHTML` 초기화.

### 내비게이션

- `search(locus)` 또는 `goto(locus)` — 유전자 조회 후 `bed_regions`에서 만든 locus 또는 변이 `chrom:pos±pad`로 이동.
- 사용자에게 “**Load / refresh IGV** 후 유전자/좌표로 점프” 안내.

---

## 5. (선택) 정적 igv-reports HTML

**데이터**: 리뷰 결과 JSON에 `igv_report_html` 같은 **상대 경로** (예 `snapshots/foo.html`).

**UI**: `fetch(orderFileUrl)` → HTML 텍스트 → (선택) DOM 파싱으로 특정 섹션 제거 → `Blob` + `createObjectURL` → `<iframe src="blob:...">`  
(직접 `file://` URL로 열기보다 **CORS/인증**이 맞는 API 경로로 fetch하는 것이 안전함.)

**큰 모달**: 동일 blob URL을 재사용해 “확대 보기” 모달에도 동일 HTML을 iframe으로 복제 가능.

---

## 6. (선택) 파이프라인 입력용 BAM CSV 브라우저

- 주문 생성 시 `sample_id,bam,bai` CSV 경로를 넣는다면, 서버에 **디렉터리 브라우징** + **CSV에서 sample_id 파싱** API를 두어 포털에서 경로를 고르게 할 수 있다 (`/api/portal/bam-csv/browse`, `sample-ids` 등).
- 이는 **IGV 뷰어와 무관**하지만 “BAM 입력” UX를 한 세트로 묶을 때 유용하다.

---

## 7. 구현 시 검증 체크리스트

- [ ] BAM/BAI가 서버에서 주문 루트 밖으로 탈출하지 않는 경로 검증
- [ ] `.bai` 인접 존재 또는 명시적 `index_rel_path`
- [ ] HEAD + GET Range (또는 동등)로 IGV.js가 대용량 BAM을 읽는지
- [ ] API 키 사용 시 IGV 트랙 `headers` 및 CORS
- [ ] ancillary BAM만 있을 때 메시지와 `params.igv_bam` 우회
- [ ] 주문 변경 시 IGV 인스턴스 dispose·메모리(blob URL revoke)

---

## 이 저장소와의 대응 관계

이 프롬프트는 다음에 대응합니다.

| 영역 | 참고 위치 (예시) |
|------|------------------|
| 백엔드 | `app/main.py`: `get_order_coverage_context`, `_list_order_bam_tracks`, `*_order_*_roots`, `head_order_file`, `download_order_file` |
| 프론트 | `portal/index.html`: `startCoverageIgvBrowser`, `coverageOrderFileUrl`, `ensureIgvScriptLoaded` 등 |

다른 프로젝트에서는 **주문 루트 정의**와 **서비스 코드별 `interpretation_genes`**만 맞추면 나머지 패턴은 그대로 이식할 수 있습니다.

대상 프로젝트 스택(예: React, 별도 백엔드)에 맞춰 프롬프트를 “컴포넌트 단위”로 쪼개려면 이 문서를 기준으로 요구사항을 나누면 됩니다.
