#!/usr/bin/env bash
# =============================================================================
#  setup_shared_refs.sh
#  공유 레퍼런스 레이어 초기화 스크립트
# =============================================================================
#
#  목적:
#    carrier-screening, sgNIPT, service-daemon 세 서비스가 공통으로 사용하는
#    대용량 레퍼런스 파일을 중앙 경로(/data/reference)에 구성하고,
#    각 서비스의 data/refs 디렉토리에 심볼릭 링크를 생성합니다.
#
#  실행 방법:
#    chmod +x setup_shared_refs.sh
#    sudo ./setup_shared_refs.sh [--ref-root /data/reference] [--dry-run]
#
#  옵션:
#    --ref-root PATH   : 중앙 레퍼런스 루트 경로 (기본값: /data/reference)
#    --cs-dir PATH     : carrier-screening 저장소 경로
#    --sgnipt-dir PATH : sgNIPT 저장소 경로
#    --dry-run         : 실제 파일 생성 없이 수행할 작업만 출력
#
#  전제 조건:
#    - 기존 레퍼런스 파일이 서버에 존재해야 합니다.
#    - BWA-MEM2 인덱스는 별도로 생성해야 합니다 (아래 안내 참고).
#    - VEP 캐시는 별도로 다운로드해야 합니다 (아래 안내 참고).
#
# =============================================================================

set -euo pipefail

# ── 기본값 설정 ───────────────────────────────────────────────────────────────
REF_ROOT="/data/reference"
CS_DIR="/home/ken/gx-exome"
SGNIPT_DIR="/home/ken/sgNIPT"
DRY_RUN=false

# 기존 레퍼런스 파일 위치 (서버 현재 상태 기준)
EXISTING_FASTA="/home/sam/dark_gene_pipeline/refs/GRCh38.fasta"
EXISTING_BWA_INDEX="/home/ken/gx-exome/data/refs/bwa_index"
EXISTING_GNOMAD="/data/reference/annotation/gnomad"
EXISTING_CLINVAR=""   # 설정 필요

# ── 인자 파싱 ─────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref-root)   REF_ROOT="$2"; shift 2 ;;
        --cs-dir)     CS_DIR="$2"; shift 2 ;;
        --sgnipt-dir) SGNIPT_DIR="$2"; shift 2 ;;
        --dry-run)    DRY_RUN=true; shift ;;
        *) echo "[WARN] 알 수 없는 옵션: $1"; shift ;;
    esac
done

# ── 헬퍼 함수 ─────────────────────────────────────────────────────────────────
log_info()  { echo -e "\033[0;32m[INFO]\033[0m  $*"; }
log_warn()  { echo -e "\033[0;33m[WARN]\033[0m  $*"; }
log_error() { echo -e "\033[0;31m[ERROR]\033[0m $*"; }
log_step()  { echo -e "\n\033[1;34m══ $* ══\033[0m"; }

run_cmd() {
    if $DRY_RUN; then
        echo "  [DRY-RUN] $*"
    else
        eval "$@"
    fi
}

make_symlink() {
    local target="$1"
    local link="$2"
    if [[ -L "$link" ]]; then
        log_warn "심볼릭 링크 이미 존재: $link → $(readlink "$link")"
    elif [[ -e "$link" ]]; then
        log_warn "파일/디렉토리 이미 존재 (링크 아님): $link — 건너뜁니다"
    else
        run_cmd ln -s "$target" "$link"
        log_info "링크 생성: $link → $target"
    fi
}

# ── 시작 배너 ─────────────────────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║         공유 레퍼런스 레이어 초기화 스크립트                    ║"
echo "║  carrier-screening + sgNIPT + service-daemon 공통 설정          ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
log_info "중앙 레퍼런스 루트 : $REF_ROOT"
log_info "carrier-screening  : $CS_DIR"
log_info "sgNIPT             : $SGNIPT_DIR"
$DRY_RUN && log_warn "DRY-RUN 모드 — 실제 파일 변경 없음"
echo ""

# =============================================================================
# STEP 1: 중앙 레퍼런스 디렉토리 구조 생성
# =============================================================================
log_step "STEP 1: 중앙 레퍼런스 디렉토리 구조 생성"

run_cmd mkdir -p "${REF_ROOT}/genomes/GRCh38"
run_cmd mkdir -p "${REF_ROOT}/genomes/GRCh38/bwa_index"
run_cmd mkdir -p "${REF_ROOT}/genomes/GRCh38/bwa_mem2_index"
run_cmd mkdir -p "${REF_ROOT}/annotation/vep_cache"
run_cmd mkdir -p "${REF_ROOT}/annotation/gnomad"
run_cmd mkdir -p "${REF_ROOT}/annotation/clinvar"
run_cmd mkdir -p "${REF_ROOT}/annotation/mane"
run_cmd mkdir -p "${REF_ROOT}/annotation/curated"
run_cmd mkdir -p "${REF_ROOT}/panels/carrier_screening"
run_cmd mkdir -p "${REF_ROOT}/panels/sgNIPT"

log_info "디렉토리 구조 생성 완료"

# =============================================================================
# STEP 2: 기존 레퍼런스 파일을 중앙 경로로 이동/링크
# =============================================================================
log_step "STEP 2: 기존 레퍼런스 파일 중앙 경로로 연결"

# GRCh38 FASTA
if [[ -f "$EXISTING_FASTA" ]]; then
    CENTRAL_FASTA="${REF_ROOT}/genomes/GRCh38/GRCh38.fasta"
    if [[ ! -f "$CENTRAL_FASTA" && ! -L "$CENTRAL_FASTA" ]]; then
        log_info "FASTA 심볼릭 링크 생성: $CENTRAL_FASTA → $EXISTING_FASTA"
        run_cmd ln -s "$EXISTING_FASTA" "$CENTRAL_FASTA"
        # .fai, .dict도 함께
        [[ -f "${EXISTING_FASTA}.fai" ]]  && run_cmd ln -s "${EXISTING_FASTA}.fai"  "${CENTRAL_FASTA}.fai"
        [[ -f "${EXISTING_FASTA%.*}.dict" ]] && run_cmd ln -s "${EXISTING_FASTA%.*}.dict" "${REF_ROOT}/genomes/GRCh38/GRCh38.dict"
    else
        log_warn "FASTA 이미 존재: $CENTRAL_FASTA"
    fi
else
    log_warn "기존 FASTA 파일 없음: $EXISTING_FASTA"
    log_warn "  → 수동으로 ${REF_ROOT}/genomes/GRCh38/GRCh38.fasta 를 설정해 주세요"
fi

# BWA-MEM 인덱스 (기존 bwa_index 디렉토리 내용 링크)
if [[ -d "$EXISTING_BWA_INDEX" ]]; then
    log_info "BWA-MEM 인덱스 파일 링크 생성"
    for f in "${EXISTING_BWA_INDEX}"/*; do
        fname=$(basename "$f")
        make_symlink "$f" "${REF_ROOT}/genomes/GRCh38/bwa_index/${fname}"
    done
else
    log_warn "기존 BWA-MEM 인덱스 없음: $EXISTING_BWA_INDEX"
fi

# gnomAD (기존 경로가 이미 /data/reference/annotation/gnomad 이면 생략)
if [[ -d "$EXISTING_GNOMAD" && "$EXISTING_GNOMAD" != "${REF_ROOT}/annotation/gnomad" ]]; then
    log_info "gnomAD 디렉토리 링크 생성"
    # 이미 중앙 경로에 있으면 그대로 사용
    if [[ ! -L "${REF_ROOT}/annotation/gnomad" ]]; then
        run_cmd rmdir "${REF_ROOT}/annotation/gnomad" 2>/dev/null || true
        run_cmd ln -s "$EXISTING_GNOMAD" "${REF_ROOT}/annotation/gnomad"
    fi
fi

log_info "STEP 2 완료"

# =============================================================================
# STEP 3: BWA-MEM2 인덱스 생성 안내
# =============================================================================
log_step "STEP 3: BWA-MEM2 인덱스 상태 확인"

BWA_MEM2_INDEX_DIR="${REF_ROOT}/genomes/GRCh38/bwa_mem2_index"
FASTA_LINK="${BWA_MEM2_INDEX_DIR}/GRCh38.fasta"
INDEX_FILE="${BWA_MEM2_INDEX_DIR}/GRCh38.fasta.bwt.2bit.64"

if [[ -f "$INDEX_FILE" ]]; then
    log_info "BWA-MEM2 인덱스 이미 존재: $INDEX_FILE"
else
    log_warn "BWA-MEM2 인덱스가 없습니다. 다음 명령으로 생성하세요:"
    echo ""
    echo "  ┌─────────────────────────────────────────────────────────────────┐"
    echo "  │  BWA-MEM2 인덱스 생성 (약 30~60분, 64GB RAM 필요)              │"
    echo "  │                                                                 │"
    echo "  │  # 1. FASTA 심볼릭 링크                                        │"
    echo "  │  ln -s ${REF_ROOT}/genomes/GRCh38/GRCh38.fasta \\"
    echo "  │        ${BWA_MEM2_INDEX_DIR}/GRCh38.fasta"
    echo "  │                                                                 │"
    echo "  │  # 2. BWA-MEM2 인덱스 생성 (Docker 사용)                       │"
    echo "  │  docker run --rm \\                                             │"
    echo "  │    -v ${BWA_MEM2_INDEX_DIR}:/refs \\                           │"
    echo "  │    quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5 \\         │"
    echo "  │    bwa-mem2 index /refs/GRCh38.fasta                           │"
    echo "  └─────────────────────────────────────────────────────────────────┘"
    echo ""
    # FASTA 링크는 미리 생성
    if [[ -f "${REF_ROOT}/genomes/GRCh38/GRCh38.fasta" || -L "${REF_ROOT}/genomes/GRCh38/GRCh38.fasta" ]]; then
        make_symlink "${REF_ROOT}/genomes/GRCh38/GRCh38.fasta" "$FASTA_LINK"
    fi
fi

# =============================================================================
# STEP 4: VEP 캐시 상태 확인
# =============================================================================
log_step "STEP 4: VEP 캐시 상태 확인"

VEP_CACHE_DIR="${REF_ROOT}/annotation/vep_cache"
VEP_SPECIES_DIR="${VEP_CACHE_DIR}/homo_sapiens"

if [[ -d "$VEP_SPECIES_DIR" ]]; then
    log_info "VEP 캐시 존재: $VEP_SPECIES_DIR"
else
    log_warn "VEP 캐시가 없습니다. 다음 명령으로 다운로드하세요 (약 15GB, 1~2시간):"
    echo ""
    echo "  ┌─────────────────────────────────────────────────────────────────┐"
    echo "  │  VEP release 111 캐시 다운로드 (carrier-screening + sgNIPT 공유) │"
    echo "  │                                                                 │"
    echo "  │  docker run --rm \\                                             │"
    echo "  │    -v ${VEP_CACHE_DIR}:/opt/vep/.vep \\                        │"
    echo "  │    ensemblorg/ensembl-vep:release_111 \\                        │"
    echo "  │    perl INSTALL.pl \\                                           │"
    echo "  │      -a cf \\                                                   │"
    echo "  │      -s homo_sapiens \\                                         │"
    echo "  │      -y GRCh38 \\                                               │"
    echo "  │      --CONVERT \\                                               │"
    echo "  │      --NO_HTSLIB                                                │"
    echo "  │                                                                 │"
    echo "  │  ※ carrier-screening과 sgNIPT가 동일한 캐시를 공유합니다.      │"
    echo "  │    한 번만 다운로드하면 됩니다.                                 │"
    echo "  └─────────────────────────────────────────────────────────────────┘"
    echo ""
fi

# =============================================================================
# STEP 5: carrier-screening 심볼릭 링크 설정
# =============================================================================
log_step "STEP 5: carrier-screening 심볼릭 링크 설정"

CS_REFS_DIR="${CS_DIR}/data/refs"
run_cmd mkdir -p "$CS_REFS_DIR"

# FASTA
make_symlink "${REF_ROOT}/genomes/GRCh38/GRCh38.fasta"      "${CS_REFS_DIR}/GRCh38.fasta"
make_symlink "${REF_ROOT}/genomes/GRCh38/GRCh38.fasta.fai"  "${CS_REFS_DIR}/GRCh38.fasta.fai"
make_symlink "${REF_ROOT}/genomes/GRCh38/GRCh38.dict"       "${CS_REFS_DIR}/GRCh38.dict"

# BWA-MEM 인덱스 (기존 호환성 유지)
make_symlink "${REF_ROOT}/genomes/GRCh38/bwa_index"          "${CS_REFS_DIR}/bwa_index"

# BWA-MEM2 인덱스 (신규)
make_symlink "${REF_ROOT}/genomes/GRCh38/bwa_mem2_index"     "${CS_REFS_DIR}/bwa_mem2_index"

# VEP 캐시
make_symlink "${REF_ROOT}/annotation/vep_cache"              "${CS_REFS_DIR}/vep_cache"

log_info "carrier-screening 링크 설정 완료"

# =============================================================================
# STEP 6: sgNIPT 심볼릭 링크 설정
# =============================================================================
log_step "STEP 6: sgNIPT 심볼릭 링크 설정"

SGNIPT_REFS_DIR="${SGNIPT_DIR}/data/refs"
run_cmd mkdir -p "$SGNIPT_REFS_DIR"

# FASTA
make_symlink "${REF_ROOT}/genomes/GRCh38/GRCh38.fasta"      "${SGNIPT_REFS_DIR}/GRCh38.fasta"
make_symlink "${REF_ROOT}/genomes/GRCh38/GRCh38.fasta.fai"  "${SGNIPT_REFS_DIR}/GRCh38.fasta.fai"
make_symlink "${REF_ROOT}/genomes/GRCh38/GRCh38.dict"       "${SGNIPT_REFS_DIR}/GRCh38.dict"

# BWA-MEM 인덱스
make_symlink "${REF_ROOT}/genomes/GRCh38/bwa_index"          "${SGNIPT_REFS_DIR}/bwa_index"

# BWA-MEM2 인덱스
make_symlink "${REF_ROOT}/genomes/GRCh38/bwa_mem2_index"     "${SGNIPT_REFS_DIR}/bwa_mem2_index"

# VEP 캐시 (carrier-screening과 동일한 캐시 공유)
make_symlink "${REF_ROOT}/annotation/vep_cache"              "${SGNIPT_REFS_DIR}/vep_cache"

log_info "sgNIPT 링크 설정 완료"

# =============================================================================
# STEP 7: Docker 컨테이너 볼륨 마운트 안내
# =============================================================================
log_step "STEP 7: Docker 볼륨 마운트 안내"

echo ""
echo "  심볼릭 링크가 Docker 컨테이너 내부에서 정상 동작하려면"
echo "  중앙 레퍼런스 경로를 컨테이너에 동일한 절대 경로로 마운트해야 합니다."
echo ""
echo "  ┌─────────────────────────────────────────────────────────────────┐"
echo "  │  carrier-screening / sgNIPT 실행 시 추가 마운트 옵션           │"
echo "  │                                                                 │"
echo "  │  docker run ... \\                                              │"
echo "  │    -v ${REF_ROOT}:${REF_ROOT}:ro \\                            │"
echo "  │    ...                                                          │"
echo "  │                                                                 │"
echo "  │  service-daemon docker-compose.yml (이미 설정됨):              │"
echo "  │    volumes:                                                     │"
echo "  │      - \${BASE_DIR:-/data}/reference:/data/reference:ro        │"
echo "  └─────────────────────────────────────────────────────────────────┘"
echo ""

# =============================================================================
# STEP 8: 서버 분리(Migration) 시나리오 안내
# =============================================================================
log_step "STEP 8: 서버 분리 시나리오 안내"

echo ""
echo "  추후 서비스를 다른 서버로 분리할 때의 절차:"
echo ""
echo "  1. 새 서버에 동일한 중앙 레퍼런스 구조 복사:"
echo "     rsync -avz --progress ${REF_ROOT}/ newserver:${REF_ROOT}/"
echo ""
echo "  2. 새 서버에서 이 스크립트 실행:"
echo "     ./setup_shared_refs.sh --ref-root ${REF_ROOT} \\"
echo "       --cs-dir /path/to/carrier-screening \\"
echo "       --sgnipt-dir /path/to/sgNIPT"
echo ""
echo "  3. 소스 코드 수정 없이 즉시 동작합니다."
echo "     (nextflow.config의 레퍼런스 경로가 /data/reference로 통일되어 있음)"
echo ""

# =============================================================================
# 완료
# =============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║  초기화 완료                                                     ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
echo "  다음 단계:"
echo "  1. BWA-MEM2 인덱스 생성 (STEP 3 안내 참고)"
echo "  2. VEP 캐시 다운로드 (STEP 4 안내 참고)"
echo "  3. carrier-screening 실행 테스트:"
echo "     nextflow run bin/main.nf --aligner bwa-mem2 --variant_caller deepvariant"
echo "  4. sgNIPT 실행 테스트:"
echo "     nextflow run bin/main.nf --aligner bwa-mem2 --skip_vep false"
echo ""
