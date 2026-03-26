#!/usr/bin/env bash
# 컨테이너 또는 호스트에서 1회 실행: snpEff 유전체 DB를 SNPEFF_DATA_DIR 에 받습니다.
# 예: docker compose exec service-daemon bash /app/scripts/snpeff-download.sh
set -euo pipefail
JAR="${SNPEFF_JAR:-/opt/tools/snpEff/snpEff.jar}"
DB="${SNPEFF_DB:-GRCh38.86}"
DD="${SNPEFF_DATA_DIR:-/data/carrier_screening_work/snpeff}"
mkdir -p "$DD"
exec java -jar "$JAR" download -v "$DB" -dataDir "$DD"
