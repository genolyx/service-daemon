# Multi-service genomics API + Test Portal (FastAPI).
# Layout aligned with nipt-daemon: non-root user, optional Docker CLI (Nextflow -with-docker), TZ.

FROM python:3.11-slim

ARG TZ=Asia/Seoul
ARG HOST_UID=1000
ARG HOST_GID=1000
ARG DOCKER_GID=999

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=${TZ} \
    PYTHONUNBUFFERED=1

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
      bash ca-certificates curl gnupg tzdata procps \
      openjdk-21-jre-headless \
      libglib2.0-0 libpango-1.0-0 libpangocairo-1.0-0 \
      libcairo2 libffi8 libfontconfig1 fonts-liberation; \
    install -m 0755 -d /etc/apt/keyrings; \
    curl -fsSL https://download.docker.com/linux/debian/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg; \
    chmod a+r /etc/apt/keyrings/docker.gpg; \
    . /etc/os-release; \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/debian ${VERSION_CODENAME} stable" > /etc/apt/sources.list.d/docker.list; \
    apt-get update; \
    apt-get install -y --no-install-recommends docker-ce-cli; \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime; echo $TZ > /etc/timezone; \
    rm -rf /var/lib/apt/lists/*

# Nextflow (컨테이너 내 실행 가능한 바이너리; 호스트 바인드 마운트는 exit 126 흔함)
RUN set -eux; \
    cd /tmp; \
    curl -fsSL https://get.nextflow.io | bash; \
    install -m 0755 nextflow /usr/local/bin/nextflow; \
    rm -f nextflow; \
    /usr/local/bin/nextflow -version

RUN groupadd -g ${HOST_GID} ken || true && \
    useradd -u ${HOST_UID} -g ${HOST_GID} -m -s /bin/bash ken || true && \
    groupadd -g ${DOCKER_GID} docker || true && \
    usermod -aG docker ken

# snpEff: 빌드 시 네트워크 없는 환경을 위해 다운로드하지 않음.
# 두 가지 방식 중 하나 선택:
#   A) 호스트에 있는 snpEff 디렉토리를 볼륨 마운트  (SNPEFF_JAR=/opt/tools/snpEff/snpEff.jar)
#   B) 컨테이너 최초 기동 후 수동 설치:
#        docker compose exec service-daemon bash scripts/snpeff-install.sh
RUN apt-get update && \
    apt-get install -y --no-install-recommends unzip && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p /opt/tools/snpEff && \
    chown -R ken:ken /opt/tools

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app/ ./app/
COPY portal/ ./portal/
COPY data/ ./data/
COPY scripts/ ./scripts/
# tools/snpEff/ — git 미추적: snpEff.jar + snpEff.config (배포 zip에서 동일 디렉터리 복사)
# JAR만 있고 config가 없으면: Cannot read config file 'snpEff.config'
COPY --chown=ken:ken tools/ /opt/tools/

RUN mkdir -p /app/logs /data/fastq /data/analysis /data/output /data/log && \
    chown -R ken:ken /app /data

EXPOSE 8000

USER ken

CMD ["/bin/sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port ${APP_PORT:-8000} --no-access-log"]
