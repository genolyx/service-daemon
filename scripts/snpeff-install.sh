#!/usr/bin/env bash
# snpEff JAR 설치 스크립트 (컨테이너 내부에서 1회 실행)
#
# 사용법:
#   docker compose exec service-daemon bash /app/scripts/snpeff-install.sh
#
# 완료 후 .env.docker 에 다음 설정 확인:
#   SNPEFF_JAR=/opt/tools/snpEff/snpEff.jar
#   SNPEFF_DB=GRCh38.86
#   SNPEFF_DATA_DIR=/data/gx-exome-work/snpeff
set -euo pipefail

INSTALL_DIR="${SNPEFF_INSTALL_DIR:-/opt/tools}"
JAR_PATH="$INSTALL_DIR/snpEff/snpEff.jar"

CFG_PATH="$INSTALL_DIR/snpEff/snpEff.config"
if [ -f "$JAR_PATH" ] && [ -f "$CFG_PATH" ]; then
    echo "[snpeff-install] Already installed: $JAR_PATH"
    java -jar "$JAR_PATH" -version 2>&1 | head -1
    exit 0
fi

echo "[snpeff-install] Downloading snpEff to $INSTALL_DIR ..."
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# 다운로드 시도 (공식 URL → GitHub 미러 순)
URLS=(
    "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"
    "https://github.com/pcingola/SnpEff/releases/latest/download/snpEff_latest_core.zip"
)

for URL in "${URLS[@]}"; do
    echo "[snpeff-install] Trying: $URL"
    if curl -fsSL --connect-timeout 30 -o snpEff.zip "$URL"; then
        echo "[snpeff-install] Download OK"
        break
    fi
    echo "[snpeff-install] Failed, trying next..."
    rm -f snpEff.zip
done

if [ ! -f snpEff.zip ]; then
    echo "[snpeff-install] ERROR: All download URLs failed."
    echo "  수동으로 snpEff_latest_core.zip 을 $INSTALL_DIR 에 복사한 후 다시 실행하세요."
    exit 1
fi

unzip -q snpEff.zip
rm -f snpEff.zip

if [ ! -f "$JAR_PATH" ]; then
    echo "[snpeff-install] ERROR: JAR not found at $JAR_PATH after unzip"
    exit 1
fi
if [ ! -f "$CFG_PATH" ]; then
    echo "[snpeff-install] ERROR: snpEff.config missing next to JAR (need full snpEff zip)"
    exit 1
fi

echo "[snpeff-install] Installed: $JAR_PATH (+ snpEff.config)"
java -jar "$JAR_PATH" -version 2>&1 | head -1

echo ""
echo "[snpeff-install] Done. DB 다운로드는 snpeff-download.sh 를 실행하세요:"
echo "  docker compose exec service-daemon bash /app/scripts/snpeff-download.sh"
