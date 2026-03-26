#!/usr/bin/env bash
# 개발용: 소스 마운트 + uvicorn --reload (docker-compose-dev.yml)
set -euo pipefail
cd "$(dirname "$0")"

if [[ ! -f .env.compose ]]; then
  echo "Missing .env.compose — copy and edit:"
  echo "  cp .env.compose.example .env.compose"
  echo "Set HOST_UID, HOST_GID, DOCKER_GID from: id -u; id -g; getent group docker"
  exit 1
fi

exec docker compose --env-file .env.compose -f docker-compose-dev.yml up --build "$@"
