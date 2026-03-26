"""
Platform API Client

플랫폼과의 모든 통신을 담당합니다.
기존 nipt-daemon의 notifier.py + aws_client.py를 일반화하여,
service_code를 포함한 API 호출을 수행합니다.
"""

import os
import json
import logging
from typing import Optional, Dict, Any, List

from .config import settings
from .models import NotificationResult, NotificationStatus, Job, OutputFile
from .auth_client import auth_request

logger = logging.getLogger(__name__)


class PlatformClient:
    """플랫폼 API 클라이언트"""

    def __init__(self):
        self.base_url = settings.platform_api_base.rstrip("/")

    @staticmethod
    def _skipped(msg: str = "Platform API disabled") -> NotificationResult:
        return NotificationResult(
            status=NotificationStatus.SKIPPED,
            message=msg,
        )

    # ─── Order 관리 ────────────────────────────────────────

    async def get_pending_orders(self, service_code: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        플랫폼에서 대기 중인 주문 목록을 가져옵니다.
        
        Args:
            service_code: 특정 서비스만 조회 (None이면 전체)
        """
        try:
            url = f"{self.base_url}/analysis/orders/pending"
            params = {}
            if service_code:
                params["service_code"] = service_code

            if not settings.platform_api_enabled:
                return []

            response = await auth_request(method="GET", url=url, params=params, timeout=30.0)
            data = response.json()

            orders = data.get("data", data.get("orders", []))
            logger.info(f"Fetched {len(orders)} pending orders" +
                        (f" for {service_code}" if service_code else ""))
            return orders

        except Exception as e:
            logger.error(f"Failed to fetch pending orders: {e}")
            return []

    async def get_order_detail(self, order_id: str) -> Optional[Dict[str, Any]]:
        """주문 상세 정보 조회"""
        try:
            if not settings.platform_api_enabled:
                return None
            url = f"{self.base_url}/analysis/order/{order_id}"
            response = await auth_request(method="GET", url=url, timeout=30.0)
            data = response.json()
            return data.get("data", data)
        except Exception as e:
            logger.error(f"Failed to get order detail for {order_id}: {e}")
            return None

    # ─── 상태 업데이트 ─────────────────────────────────────

    async def update_order_status(
        self,
        order_id: str,
        service_code: str,
        status: str,
        progress: int = 0,
        message: str = ""
    ) -> NotificationResult:
        """
        주문 상태를 플랫폼에 업데이트합니다.
        
        Args:
            order_id: 주문 ID
            service_code: 서비스 코드
            status: 상태 문자열
            progress: 진행률 (0-100)
            message: 상태 메시지
        """
        try:
            if not settings.platform_api_enabled:
                return self._skipped()
            url = f"{self.base_url}/analysis/order/{order_id}/status"
            payload = {
                "service_code": service_code,
                "status": status,
                "progress": progress,
                "message": message
            }

            response = await auth_request(
                method="PUT",
                url=url,
                json=payload,
                timeout=30.0
            )

            if response.status_code == 200:
                logger.info(
                    f"[{service_code}] Updated order {order_id}: "
                    f"{status} ({progress}%) - {message}"
                )
                return NotificationResult(
                    status=NotificationStatus.SUCCESS,
                    message="Status updated",
                    response_code=response.status_code
                )
            else:
                logger.error(
                    f"[{service_code}] Status update failed for {order_id}: "
                    f"{response.status_code}"
                )
                return NotificationResult(
                    status=NotificationStatus.FAILED,
                    message=f"Status update failed: {response.status_code}",
                    response_code=response.status_code
                )

        except Exception as e:
            logger.error(f"[{service_code}] Failed to update status for {order_id}: {e}")
            return NotificationResult(
                status=NotificationStatus.FAILED,
                message=str(e)
            )

    # ─── 결과 알림 ─────────────────────────────────────────

    async def notify_analysis_result(
        self,
        order_id: str,
        service_code: str,
        success: bool,
        log: str = ""
    ) -> NotificationResult:
        """분석 결과 알림"""
        try:
            if not settings.platform_api_enabled:
                return self._skipped()
            url = f"{self.base_url}/analysis/order/{order_id}/result"
            payload = {
                "service_code": service_code,
                "success": success,
                "log": log
            }

            response = await auth_request(
                method="POST",
                url=url,
                json=payload,
                timeout=30.0
            )

            if response.status_code == 200:
                status_str = "SUCCESS" if success else "FAILED"
                logger.info(f"[{service_code}] Notified result for {order_id}: {status_str}")
                return NotificationResult(
                    status=NotificationStatus.SUCCESS,
                    message=f"Result notification sent: {status_str}",
                    response_code=response.status_code
                )
            else:
                logger.error(
                    f"[{service_code}] Result notification failed for {order_id}: "
                    f"{response.status_code}"
                )
                return NotificationResult(
                    status=NotificationStatus.FAILED,
                    message=f"Notification failed: {response.status_code}",
                    response_code=response.status_code
                )

        except Exception as e:
            logger.error(f"[{service_code}] Failed to notify result for {order_id}: {e}")
            return NotificationResult(
                status=NotificationStatus.FAILED,
                message=str(e)
            )

    async def notify_analysis_failed(
        self,
        order_id: str,
        service_code: str,
        failed_reason: str
    ) -> NotificationResult:
        """분석 실패 알림"""
        try:
            if not settings.platform_api_enabled:
                return self._skipped()
            url = f"{self.base_url}/analysis/order/{order_id}/failed"
            payload = {
                "service_code": service_code,
                "failed_reason": failed_reason
            }

            response = await auth_request(
                method="POST",
                url=url,
                json=payload,
                timeout=30.0
            )

            if response.status_code == 200:
                logger.info(f"[{service_code}] Notified failure for {order_id}")
                return NotificationResult(
                    status=NotificationStatus.SUCCESS,
                    message="Failure notification sent",
                    response_code=response.status_code
                )
            else:
                return NotificationResult(
                    status=NotificationStatus.FAILED,
                    message=f"Failure notification failed: {response.status_code}",
                    response_code=response.status_code
                )

        except Exception as e:
            logger.error(f"[{service_code}] Failed to notify failure for {order_id}: {e}")
            return NotificationResult(
                status=NotificationStatus.FAILED,
                message=str(e)
            )

    # ─── 파일 업로드 ───────────────────────────────────────

    async def upload_output_file(
        self,
        order_id: str,
        service_code: str,
        output_file: OutputFile
    ) -> NotificationResult:
        """결과 파일을 플랫폼에 업로드합니다."""
        try:
            if not settings.platform_api_enabled:
                return self._skipped()
            if not os.path.exists(output_file.file_path):
                error_msg = f"File not found: {output_file.file_path}"
                logger.error(f"[{service_code}] {error_msg}")
                return NotificationResult(
                    status=NotificationStatus.NOT_FOUND,
                    message=error_msg
                )

            file_size = os.path.getsize(output_file.file_path)
            file_name = output_file.file_name or os.path.basename(output_file.file_path)

            logger.info(
                f"[{service_code}] Uploading {output_file.file_type} for {order_id}: "
                f"{file_name} ({file_size} bytes)"
            )

            url = f"{self.base_url}/analysis/order/{order_id}/file"

            with open(output_file.file_path, "rb") as f:
                files = {"file": (file_name, f, output_file.content_type)}
                data = {
                    "service_code": service_code,
                    "file_type": output_file.file_type
                }

                response = await auth_request(
                    method="POST",
                    url=url,
                    files=files,
                    data=data,
                    timeout=300.0
                )

            if response.status_code == 200:
                logger.info(
                    f"[{service_code}] Successfully uploaded {output_file.file_type} "
                    f"for {order_id}"
                )
                return NotificationResult(
                    status=NotificationStatus.SUCCESS,
                    message=f"File uploaded: {file_name}",
                    response_code=response.status_code
                )
            else:
                logger.error(
                    f"[{service_code}] File upload failed for {order_id}: "
                    f"{response.status_code} - {response.text}"
                )
                return NotificationResult(
                    status=NotificationStatus.FAILED,
                    message=f"Upload failed: {response.status_code}",
                    response_code=response.status_code
                )

        except Exception as e:
            logger.error(f"[{service_code}] Error uploading file for {order_id}: {e}")
            return NotificationResult(
                status=NotificationStatus.FAILED,
                message=str(e)
            )

    async def upload_pdf_report(
        self,
        order_id: str,
        service_code: str,
        pdf_path: str,
        signed: bool = False
    ) -> NotificationResult:
        """PDF 리포트 업로드"""
        try:
            if not settings.platform_api_enabled:
                return self._skipped()
            if not os.path.exists(pdf_path):
                return NotificationResult(
                    status=NotificationStatus.NOT_FOUND,
                    message=f"PDF file not found: {pdf_path}"
                )

            file_size = os.path.getsize(pdf_path)
            kind = "SIGNED" if signed else "NORMAL"
            logger.info(
                f"[{service_code}] Uploading {kind} PDF report for {order_id}: "
                f"{file_size} bytes ({pdf_path})"
            )

            if signed:
                url = f"{self.base_url}/analysis/order/{order_id}/signed-pdf-result"
            else:
                url = f"{self.base_url}/analysis/order/{order_id}/pdf-result"

            with open(pdf_path, "rb") as f:
                files = {"file": f}
                data = {"service_code": service_code}

                response = await auth_request(
                    method="POST",
                    url=url,
                    files=files,
                    data=data,
                    timeout=120.0
                )

            if response.status_code == 200:
                logger.info(
                    f"[{service_code}] Successfully uploaded {kind} PDF report for {order_id}"
                )
                return NotificationResult(
                    status=NotificationStatus.SUCCESS,
                    message="PDF report uploaded successfully",
                    response_code=response.status_code
                )
            else:
                logger.error(
                    f"[{service_code}] PDF upload failed for {order_id}: "
                    f"{response.status_code}"
                )
                return NotificationResult(
                    status=NotificationStatus.FAILED,
                    message=f"PDF upload failed: {response.status_code}",
                    response_code=response.status_code
                )

        except Exception as e:
            logger.error(f"[{service_code}] Failed to upload PDF for {order_id}: {e}")
            return NotificationResult(
                status=NotificationStatus.FAILED,
                message=str(e)
            )

    # ─── 일괄 업로드 ───────────────────────────────────────

    async def upload_all_outputs(
        self,
        order_id: str,
        service_code: str,
        output_files: List[OutputFile]
    ) -> Dict[str, NotificationResult]:
        """모든 결과 파일을 일괄 업로드합니다."""
        if not settings.platform_api_enabled:
            sk = self._skipped()
            return {of.file_type: sk for of in output_files}
        results = {}
        for output_file in output_files:
            result = await self.upload_output_file(order_id, service_code, output_file)
            results[output_file.file_type] = result

            if result.status == NotificationStatus.FAILED:
                logger.warning(
                    f"[{service_code}] Upload failed for {output_file.file_type}, "
                    f"continuing with remaining files..."
                )

        success_count = sum(
            1 for r in results.values() if r.status == NotificationStatus.SUCCESS
        )
        logger.info(
            f"[{service_code}] Upload complete for {order_id}: "
            f"{success_count}/{len(output_files)} files uploaded"
        )
        return results


# 전역 인스턴스
_platform_client = PlatformClient()


def get_platform_client() -> PlatformClient:
    """전역 PlatformClient 인스턴스 반환"""
    return _platform_client
