#!/bin/bash
# Service Daemon 연동 테스트 실행 스크립트
# sgNIPT + carrier_screening + Test Portal
#
# 사용:
#   ./run-integration.sh              # 포그라운드 (Ctrl+C 종료)
#   ./run-integration.sh -b           # 백그라운드 (nohup + logs/integration.log)
#   INTEGRATION_LOG=/tmp/sd.log ./run-integration.sh --background

set -e
cd "$(dirname "$0")"

BG=0
filtered=()
for a in "$@"; do
  case "$a" in
    -b|--background) BG=1 ;;
    *) filtered+=("$a") ;;
  esac
done

if [ "$BG" = 1 ]; then
  mkdir -p logs
  LOGFILE="${INTEGRATION_LOG:-$PWD/logs/integration.log}"
  nohup "$0" "${filtered[@]}" >>"$LOGFILE" 2>&1 &
  pid=$!
  echo "$pid" > .integration.pid
  echo "Started PID $pid (nohup)"
  echo "  Log: $LOGFILE   → tail -f \"$LOGFILE\""
  echo "  Stop: kill \$(cat \"$PWD/.integration.pid\")"
  exit 0
fi

set -- "${filtered[@]}"

# .env 설정 (없으면 .env.integration 복사)
if [ ! -f .env ]; then
  if [ -f .env.integration ]; then
    echo "Creating .env from .env.integration..."
    cp .env.integration .env
  else
    echo "Error: .env.integration not found"
    exit 1
  fi
fi
# .env 로드 (APP_PORT 등)
set -a
[ -f .env ] && . .env
set +a

# 가상환경 확인
if [ ! -d .venv ]; then
  echo "Creating virtual environment..."
  python3 -m venv .venv
  .venv/bin/pip install -r requirements.txt
fi

# nipt-daemon: 8000, 8001 / service-daemon: 8002, 8003...
PORT=${APP_PORT:-8003}
for p in 8003 8002 8080 8888; do
  if ! ss -tlnp 2>/dev/null | grep -q ":$p "; then
    PORT=$p
    break
  fi
done

echo "=============================================="
echo "  Service Daemon - Integration Test"
echo "  Port: $PORT"
echo "  Portal: http://localhost:$PORT/portal"
echo "  Health: http://localhost:$PORT/health"
echo "=============================================="

.venv/bin/python -m uvicorn app.main:app --host 0.0.0.0 --port "$PORT" --no-access-log
