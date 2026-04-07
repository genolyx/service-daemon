# Service Daemon
This project will contain the generalized service daemon.

## 연동 테스트 (sgNIPT + carrier_screening + Portal)

1. **환경 설정**
   ```bash
   cp .env.integration .env
   # .env에서 CARRIER_SCREENING_PIPELINE_DIR, SGNIPT_PIPELINE_DIR 등 경로 확인
   # Portal에서 Whole Exome / Health Screening 주문을 쓰려면 ENABLED_SERVICES에
   # whole_exome, health_screening가 포함되어야 합니다 (예시 .env.integration 참고).
   ```

2. **실행**
   ```bash
   ./run-integration.sh
   ```
   (nipt-daemon: 8000/8001, service-daemon: 8002 또는 8003)

3. **접속**
   - Portal: http://localhost:8003/portal (실행 포트에 맞게)
   - Health: http://localhost:8003/health
   - API Docs: http://localhost:8003/docs

4. **Portal에서 Daemon URL 변경**: Settings 버튼 → Daemon URL 입력 (다른 서버/포트인 경우)
