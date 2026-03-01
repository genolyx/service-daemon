"""
Pipeline Runner

큐에서 작업을 가져와 서비스 플러그인을 통해 파이프라인을 실행합니다.
기존 nipt-daemon의 runner.py를 일반화하여 플러그인 기반으로 동작합니다.
"""

import os
import asyncio
import logging
import signal
from datetime import datetime
from typing import Optional

from .config import settings
from .models import Job, OrderStatus, NotificationStatus
from .queue_manager import get_queue_manager
from .platform_client import get_platform_client
from .services import get_plugin

logger = logging.getLogger(__name__)


class PipelineRunner:
    """
    범용 파이프라인 실행기.
    
    큐에서 작업을 가져와 해당 서비스의 플러그인을 통해
    파이프라인을 실행하고 결과를 처리합니다.
    """

    def __init__(self):
        self._queue_manager = get_queue_manager()
        self._platform_client = get_platform_client()
        self._shutdown_event = asyncio.Event()
        self._active_processes: dict = {}  # order_id -> asyncio.subprocess.Process

    async def start(self):
        """워커 루프 시작"""
        logger.info(
            f"PipelineRunner started "
            f"(max_concurrent={self._queue_manager.max_concurrent})"
        )

        workers = [
            asyncio.create_task(self._worker(i))
            for i in range(self._queue_manager.max_concurrent)
        ]

        # Graceful shutdown 시그널 등록
        loop = asyncio.get_event_loop()
        for sig in (signal.SIGTERM, signal.SIGINT):
            try:
                loop.add_signal_handler(sig, self._handle_shutdown)
            except NotImplementedError:
                pass  # Windows에서는 지원되지 않음

        await self._shutdown_event.wait()

        # 워커 정리
        for worker in workers:
            worker.cancel()
        await asyncio.gather(*workers, return_exceptions=True)

        logger.info("PipelineRunner stopped")

    def _handle_shutdown(self):
        """Graceful shutdown 처리"""
        logger.info("Shutdown signal received, stopping workers...")
        self._shutdown_event.set()

    async def _worker(self, worker_id: int):
        """개별 워커 루프"""
        logger.info(f"Worker-{worker_id} started")

        while not self._shutdown_event.is_set():
            try:
                # 큐에서 작업 가져오기 (타임아웃 포함)
                try:
                    job = await asyncio.wait_for(
                        self._queue_manager.dequeue(),
                        timeout=settings.queue_poll_interval
                    )
                except asyncio.TimeoutError:
                    continue

                # 실행 슬롯 획득
                await self._queue_manager.acquire_slot()

                try:
                    await self._execute_job(job, worker_id)
                finally:
                    self._queue_manager.release_slot()

            except asyncio.CancelledError:
                break
            except Exception as e:
                logger.error(f"Worker-{worker_id} unexpected error: {e}", exc_info=True)
                await asyncio.sleep(5)

        logger.info(f"Worker-{worker_id} stopped")

    async def _execute_job(self, job: Job, worker_id: int):
        """
        단일 작업 실행 (전체 라이프사이클).
        
        1. 플러그인 조회
        2. 입력 준비
        3. 파이프라인 실행
        4. 완료 확인
        5. 결과 후처리
        6. 파일 업로드
        7. 완료 알림
        """
        plugin = get_plugin(job.service_code)
        if not plugin:
            error_msg = f"No plugin found for service_code: {job.service_code}"
            logger.error(error_msg)
            await self._queue_manager.mark_failed(job, error_msg)
            await self._platform_client.notify_analysis_failed(
                job.order_id, job.service_code, error_msg
            )
            return

        logger.info(
            f"Worker-{worker_id} executing job {job.order_id} "
            f"[{plugin.display_name}] (sample: {job.sample_name})"
        )

        await self._queue_manager.mark_running(job)
        await plugin.on_job_start(job)

        try:
            # ── Step 1: 입력 준비 ──
            await self._update_status(job, OrderStatus.DOWNLOADING, 5, "Preparing inputs...")
            input_ok = await plugin.prepare_inputs(job)
            if not input_ok:
                raise RuntimeError("Failed to prepare input files")

            # ── Step 2: 파이프라인 실행 ──
            await self._update_status(job, OrderStatus.RUNNING, 10, "Starting pipeline...")
            command = await plugin.get_pipeline_command(job)
            logger.info(f"[{job.service_code}] Pipeline command: {command}")

            exit_code = await self._run_pipeline(job, command)

            if exit_code != 0:
                raise RuntimeError(
                    f"Pipeline exited with code {exit_code}"
                )

            # ── Step 3: 완료 확인 ──
            await self._update_status(job, OrderStatus.RUNNING, 85, "Verifying results...")
            completion_ok = await plugin.check_completion(job)
            if not completion_ok:
                raise RuntimeError("Pipeline completion check failed")

            # ── Step 4: 결과 후처리 ──
            await self._update_status(job, OrderStatus.PROCESSING, 90, "Processing results...")
            process_ok = await plugin.process_results(job)
            if not process_ok:
                raise RuntimeError("Result processing failed")

            # ── Step 5: 파일 업로드 ──
            await self._update_status(job, OrderStatus.UPLOADING, 95, "Uploading results...")
            output_files = await plugin.get_output_files(job)

            if output_files:
                upload_results = await self._platform_client.upload_all_outputs(
                    job.order_id, job.service_code, output_files
                )
                
                # 업로드 실패 확인
                failed_uploads = [
                    ft for ft, r in upload_results.items()
                    if r.status == NotificationStatus.FAILED
                ]
                if failed_uploads:
                    logger.warning(
                        f"[{job.service_code}] Some uploads failed: {failed_uploads}"
                    )

            # ── Step 6: 완료 알림 ──
            await self._update_status(job, OrderStatus.COMPLETED, 100, "Analysis completed")
            await self._platform_client.notify_analysis_result(
                job.order_id, job.service_code, success=True
            )
            await self._queue_manager.mark_completed(job)
            await plugin.on_job_complete(job)

        except Exception as e:
            error_msg = str(e)
            logger.error(
                f"[{job.service_code}] Job {job.order_id} failed: {error_msg}",
                exc_info=True
            )
            await self._queue_manager.mark_failed(job, error_msg)
            await self._platform_client.notify_analysis_failed(
                job.order_id, job.service_code, error_msg
            )
            await plugin.on_job_failed(job, error_msg)

        finally:
            # 프로세스 참조 정리
            self._active_processes.pop(job.order_id, None)
            
            # 플러그인 정리
            try:
                await plugin.cleanup(job)
            except Exception as e:
                logger.warning(f"[{job.service_code}] Cleanup error: {e}")

    async def _run_pipeline(self, job: Job, command: str) -> int:
        """
        파이프라인 프로세스를 실행하고 모니터링합니다.
        
        Returns:
            프로세스 exit code
        """
        # 로그 디렉토리 생성
        log_dir = job.log_dir or os.path.join(
            settings.log_base_dir, job.service_code, job.work_dir, job.sample_name
        )
        os.makedirs(log_dir, exist_ok=True)
        
        log_file = os.path.join(log_dir, "pipeline.log")
        
        logger.info(f"[{job.service_code}] Running pipeline, log: {log_file}")

        with open(log_file, "w") as log_f:
            process = await asyncio.create_subprocess_shell(
                command,
                stdout=log_f,
                stderr=asyncio.subprocess.STDOUT,
                cwd=job.analysis_dir or settings.analysis_base_dir
            )

        job.pid = process.pid
        self._active_processes[job.order_id] = process

        logger.info(f"[{job.service_code}] Pipeline PID: {process.pid}")

        # 진행률 모니터링 태스크
        monitor_task = asyncio.create_task(
            self._monitor_progress(job, log_file)
        )

        # 프로세스 완료 대기
        exit_code = await process.wait()
        job.exit_code = exit_code

        # 모니터링 태스크 정리
        monitor_task.cancel()
        try:
            await monitor_task
        except asyncio.CancelledError:
            pass

        logger.info(
            f"[{job.service_code}] Pipeline finished with exit code: {exit_code}"
        )
        return exit_code

    async def _monitor_progress(self, job: Job, log_file: str):
        """
        파이프라인 로그를 모니터링하여 진행률을 업데이트합니다.
        """
        plugin = get_plugin(job.service_code)
        if not plugin:
            return

        progress_stages = plugin.get_progress_stages()
        if not progress_stages:
            return

        last_progress = 10  # 시작 진행률

        try:
            while True:
                await asyncio.sleep(30)  # 30초마다 체크

                if not os.path.exists(log_file):
                    continue

                try:
                    current_progress = self._parse_progress(log_file, progress_stages)
                    if current_progress > last_progress:
                        last_progress = current_progress
                        await self._update_status(
                            job, OrderStatus.RUNNING, current_progress,
                            f"Pipeline running ({current_progress}%)"
                        )
                except Exception as e:
                    logger.debug(f"Progress parsing error: {e}")

        except asyncio.CancelledError:
            pass

    def _parse_progress(self, log_file: str, stages: dict) -> int:
        """로그 파일에서 진행률을 파싱합니다."""
        progress = 0
        try:
            with open(log_file, "r") as f:
                for line in f:
                    for stage_name, stage_pct in stages.items():
                        if stage_name in line and stage_pct > progress:
                            progress = stage_pct
        except Exception:
            pass
        return progress

    async def _update_status(
        self, job: Job, status: OrderStatus, progress: int, message: str
    ):
        """작업 상태를 내부 + 플랫폼에 동시 업데이트"""
        job.update_status(status, progress, message)
        
        # 플랫폼에 비동기 업데이트 (실패해도 계속 진행)
        try:
            await self._platform_client.update_order_status(
                job.order_id, job.service_code,
                status.value, progress, message
            )
        except Exception as e:
            logger.warning(
                f"[{job.service_code}] Failed to update platform status: {e}"
            )

    async def cancel_job(self, order_id: str) -> bool:
        """실행 중인 작업을 취소합니다."""
        process = self._active_processes.get(order_id)
        if process and process.returncode is None:
            try:
                process.send_signal(signal.SIGTERM)
                logger.info(f"Sent SIGTERM to job {order_id} (PID: {process.pid})")
                
                # 5초 후에도 종료되지 않으면 SIGKILL
                try:
                    await asyncio.wait_for(process.wait(), timeout=5.0)
                except asyncio.TimeoutError:
                    process.kill()
                    logger.warning(f"Force killed job {order_id} (PID: {process.pid})")
                
                return True
            except Exception as e:
                logger.error(f"Failed to cancel job {order_id}: {e}")
                return False
        
        logger.warning(f"No active process found for job {order_id}")
        return False


# 전역 인스턴스
_runner: Optional[PipelineRunner] = None


def get_runner() -> PipelineRunner:
    """전역 PipelineRunner 인스턴스 반환"""
    global _runner
    if _runner is None:
        _runner = PipelineRunner()
    return _runner
