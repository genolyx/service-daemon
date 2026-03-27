#!/usr/bin/env bash
# =============================================================================
#  export_service.sh
#  서비스 독립 배포 스크립트 (서버 분리 시 사용)
# =============================================================================
#
#  목적:
#    현재 서버에서 특정 서비스(carrier-screening 또는 sgNIPT)를
#    독립적으로 다른 서버에 배포할 때 필요한 파일 패키지를 생성합니다.
#
#  사용 방법:
#    chmod +x export_service.sh
#    ./export_service.sh --service carrier-screening --dest /tmp/export
#    ./export_service.sh --service sgnipt --dest /tmp/export
#    ./export_service.sh --service all --dest /tmp/export
#
#  옵션:
#    --service NAME  : 배포할 서비스 (carrier-screening | sgnipt | all)
#    --dest PATH     : 내보낼 대상 디렉토리
#    --ref-root PATH : 중앙 레퍼런스 루트 (기본값: /data/reference)
#    --skip-refs     : 레퍼런스 파일 제외 (소스 코드만 패키징)
#
# =============================================================================

set -euo pipefail

# ── 기본값 ────────────────────────────────────────────────────────────────────
SERVICE="all"
DEST_DIR="/tmp/genolyx_export"
REF_ROOT="/data/reference"
SKIP_REFS=false

# ── 인자 파싱 ─────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --service)    SERVICE="$2"; shift 2 ;;
        --dest)       DEST_DIR="$2"; shift 2 ;;
        --ref-root)   REF_ROOT="$2"; shift 2 ;;
        --skip-refs)  SKIP_REFS=true; shift ;;
        *) echo "[WARN] 알 수 없는 옵션: $1"; shift ;;
    esac
done

log_info()  { echo -e "\033[0;32m[INFO]\033[0m  $*"; }
log_warn()  { echo -e "\033[0;33m[WARN]\033[0m  $*"; }
log_step()  { echo -e "\n\033[1;34m══ $* ══\033[0m"; }

echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║         서비스 독립 배포 패키지 생성                            ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
log_info "서비스   : $SERVICE"
log_info "대상 경로: $DEST_DIR"
log_info "레퍼런스 : $REF_ROOT"
$SKIP_REFS && log_warn "레퍼런스 파일 제외 모드"

mkdir -p "$DEST_DIR"

# =============================================================================
# 공통: 레퍼런스 파일 목록 생성 (rsync manifest)
# =============================================================================
log_step "레퍼런스 파일 목록 생성"

cat > "${DEST_DIR}/ref_manifest.txt" << EOF
# =============================================================================
#  레퍼런스 파일 목록 (새 서버로 rsync 시 사용)
# =============================================================================
#
#  새 서버로 복사 명령:
#    rsync -avz --progress --files-from=ref_manifest.txt \
#          sourceserver:/ /
#
# =============================================================================

# ── 필수: GRCh38 FASTA ──────────────────────────────────────────────────────
${REF_ROOT}/genomes/GRCh38/GRCh38.fasta
${REF_ROOT}/genomes/GRCh38/GRCh38.fasta.fai
${REF_ROOT}/genomes/GRCh38/GRCh38.dict

# ── 필수: BWA-MEM 인덱스 ────────────────────────────────────────────────────
${REF_ROOT}/genomes/GRCh38/bwa_index/

# ── 권장: BWA-MEM2 인덱스 ───────────────────────────────────────────────────
${REF_ROOT}/genomes/GRCh38/bwa_mem2_index/

# ── 권장: VEP 캐시 (carrier-screening + sgNIPT 공유, ~15GB) ─────────────────
${REF_ROOT}/annotation/vep_cache/

# ── 선택: gnomAD v3.1.2 genomes (VEP 사용 시 불필요) ───────────────────────
# ${REF_ROOT}/annotation/gnomad/

# ── 필수: ClinVar ────────────────────────────────────────────────────────────
${REF_ROOT}/annotation/clinvar/

# ── 필수: MANE RefSeq ────────────────────────────────────────────────────────
${REF_ROOT}/annotation/mane/

# ── 선택: Curated Variants DB ────────────────────────────────────────────────
${REF_ROOT}/annotation/curated/
EOF

log_info "레퍼런스 목록 생성: ${DEST_DIR}/ref_manifest.txt"

# =============================================================================
# 새 서버 초기화 스크립트 생성
# =============================================================================
log_step "새 서버 초기화 스크립트 생성"

cat > "${DEST_DIR}/setup_new_server.sh" << 'SETUP_EOF'
#!/usr/bin/env bash
# =============================================================================
#  새 서버 초기화 스크립트 (자동 생성됨)
# =============================================================================
set -euo pipefail

REF_ROOT="${1:-/data/reference}"
SERVICE="${2:-all}"

echo "[INFO] 새 서버 초기화 시작..."
echo "[INFO] 레퍼런스 루트: $REF_ROOT"
echo "[INFO] 서비스: $SERVICE"

# 디렉토리 구조 생성
mkdir -p "${REF_ROOT}/genomes/GRCh38/bwa_index"
mkdir -p "${REF_ROOT}/genomes/GRCh38/bwa_mem2_index"
mkdir -p "${REF_ROOT}/annotation/vep_cache"
mkdir -p "${REF_ROOT}/annotation/gnomad"
mkdir -p "${REF_ROOT}/annotation/clinvar"
mkdir -p "${REF_ROOT}/annotation/mane"
mkdir -p "${REF_ROOT}/annotation/curated"

echo "[INFO] 디렉토리 구조 생성 완료"
echo ""
echo "다음 단계:"
echo "1. 레퍼런스 파일 복사:"
echo "   rsync -avz --progress sourceserver:${REF_ROOT}/ ${REF_ROOT}/"
echo ""
echo "2. 서비스 심볼릭 링크 설정:"
echo "   ./setup_shared_refs.sh --ref-root ${REF_ROOT}"
echo ""
echo "3. BWA-MEM2 인덱스가 없으면 생성:"
echo "   docker run --rm \\"
echo "     -v ${REF_ROOT}/genomes/GRCh38/bwa_mem2_index:/refs \\"
echo "     quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5 \\"
echo "     bwa-mem2 index /refs/GRCh38.fasta"
SETUP_EOF

chmod +x "${DEST_DIR}/setup_new_server.sh"
log_info "새 서버 초기화 스크립트 생성: ${DEST_DIR}/setup_new_server.sh"

# =============================================================================
# service-daemon .env 템플릿 생성 (독립 배포용)
# =============================================================================
log_step "service-daemon .env 템플릿 생성"

cat > "${DEST_DIR}/.env.standalone" << EOF
# =============================================================================
#  service-daemon 독립 배포 환경 설정 (자동 생성됨)
#  서비스: ${SERVICE}
# =============================================================================

APP_NAME=service-daemon
APP_ENV=prod
APP_PORT=8000

# ── 활성화할 서비스 ──────────────────────────────────────────────────────────
EOF

if [[ "$SERVICE" == "carrier-screening" ]]; then
    echo "ENABLED_SERVICES=carrier_screening" >> "${DEST_DIR}/.env.standalone"
elif [[ "$SERVICE" == "sgnipt" ]]; then
    echo "ENABLED_SERVICES=sgnipt" >> "${DEST_DIR}/.env.standalone"
else
    echo "ENABLED_SERVICES=carrier_screening,sgnipt" >> "${DEST_DIR}/.env.standalone"
fi

cat >> "${DEST_DIR}/.env.standalone" << EOF

# ── 기본 디렉토리 ────────────────────────────────────────────────────────────
BASE_DIR=/data
FASTQ_BASE_DIR=/data/fastq
ANALYSIS_BASE_DIR=/data/analysis
OUTPUT_BASE_DIR=/data/output
LOG_BASE_DIR=/data/log

# ── 파이프라인 경로 ──────────────────────────────────────────────────────────
CARRIER_SCREENING_PIPELINE_DIR=/opt/pipelines/carrier-screening
SGNIPT_PIPELINE_DIR=/opt/pipelines/sgnipt

# ── 공유 레퍼런스 (중앙 경로) ────────────────────────────────────────────────
REF_FASTA=${REF_ROOT}/genomes/GRCh38/GRCh38.fasta
REF_FAI=${REF_ROOT}/genomes/GRCh38/GRCh38.fasta.fai
REF_DICT=${REF_ROOT}/genomes/GRCh38/GRCh38.dict
REF_BWA_INDICES=${REF_ROOT}/genomes/GRCh38/bwa_index
REF_BWA_MEM2_INDICES=${REF_ROOT}/genomes/GRCh38/bwa_mem2_index

# ── Annotation ───────────────────────────────────────────────────────────────
CLINVAR_VCF=${REF_ROOT}/annotation/clinvar/clinvar_latest.vcf.gz
GNOMAD_DIR=${REF_ROOT}/annotation/gnomad
MANE_GFF=${REF_ROOT}/annotation/mane/MANE.GRCh38.v1.3.refseq_genomic.gff.gz
CURATED_VARIANTS_DB=${REF_ROOT}/annotation/curated/curated_variants.sqlite

# VEP 사용 시 (skip_vep=false): gnomAD 로컬 파일 불필요
# VEP 미사용 시 (skip_vep=true): gnomAD 로컬 파일 필요
GNOMAD_GENOMES_GLOB=gnomad.genomes.v3.1.2.sites.*.vcf.bgz
EOF

log_info ".env 템플릿 생성: ${DEST_DIR}/.env.standalone"

# =============================================================================
# 완료 요약
# =============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║  패키지 생성 완료                                                ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
echo "  생성된 파일:"
echo "    ${DEST_DIR}/ref_manifest.txt      ← 레퍼런스 파일 목록 (rsync용)"
echo "    ${DEST_DIR}/setup_new_server.sh   ← 새 서버 초기화 스크립트"
echo "    ${DEST_DIR}/.env.standalone       ← service-daemon 환경 설정 템플릿"
echo ""
echo "  새 서버 배포 절차:"
echo "  1. 새 서버에서: bash setup_new_server.sh"
echo "  2. 레퍼런스 복사: rsync -avz --files-from=ref_manifest.txt \\"
echo "       sourceserver:/ /"
echo "  3. 소스 코드 배포: git clone https://github.com/genolyx/carrier-screening"
echo "  4. .env.standalone → .env 복사 후 수정"
echo "  5. docker compose up -d"
echo ""
