"""
Queue Manager

작업 큐와 동시 실행 제어를 관리합니다.
기존 nipt-daemon의 queue_manager.py를 일반화하여 멀티서비스를 지원합니다.
"""

import asyncio
import logging
from typing import Dict, List, Optional, Set
from collections import defaultdict

from .config import settings
from .datetime_kst import now_kst_iso
from .models import Job, OrderStatus, QueueSummary
from .order_store import (
    OrderStore,
    ACTIVE_BEFORE_RESTART,
    ingest_result_json_from_disk,
)
from .order_cleanup import delete_run_artifacts
from .services.carrier_screening.layout_norm import apply_carrier_layout_directories
from .services.sgnipt import apply_sgnipt_layout_directories

logger = logging.getLogger(__name__)


class QueueManager:
    """
    멀티서비스 작업 큐 관리자.

    단일 큐에서 모든 서비스의 작업을 관리하며,
    Semaphore를 통해 동시 실행 수를 제어합니다.
    """

    def __init__(self, max_concurrent: int = None, store: Optional[OrderStore] = None):
        self._max_concurrent = max_concurrent or settings.max_concurrent_jobs
        self._queue: asyncio.Queue = asyncio.Queue()
        self._semaphore = asyncio.Semaphore(self._max_concurrent)
        self._store = store

        # 작업 상태 추적
        self._saved_jobs: Dict[str, Job] = {}  # order_id -> Job (SAVED, 파이프라인 미시작)
        self._jobs: Dict[str, Job] = {}  # order_id -> Job (QUEUED, 큐 대기)
        self._running_jobs: Dict[str, Job] = {}  # order_id -> Job (실행 중)
        self._completed_jobs: Dict[str, Job] = {}  # order_id -> Job (완료/실패/취소)
        self._cancel_requested: Set[str] = set()  # QUEUED 작업 dequeue 시 버림

        # 서비스별 통계
        self._stats: Dict[str, Dict[str, int]] = defaultdict(
            lambda: {"queued": 0, "running": 0, "completed": 0, "failed": 0}
        )

        self._lock = asyncio.Lock()
        if self._store:
            self._restore_from_store()
        logger.info(
            f"QueueManager initialized (max_concurrent={self._max_concurrent}, "
            f"store={'on' if self._store else 'off'})"
        )

    @property
    def store(self) -> Optional[OrderStore]:
        return self._store

    def _restore_from_store(self) -> None:
        assert self._store is not None
        jobs = self._store.fetch_all_jobs()
        if not jobs:
            logger.info("No orders to restore from SQLite")
            return
        queued_list: List[Job] = []
        interrupted = 0
        for job in jobs:
            st = job.status
            if (
                job.service_code == "carrier_screening"
                and st not in (OrderStatus.COMPLETED, OrderStatus.REPORT_READY)
                and apply_carrier_layout_directories(job)
            ):
                self._store.upsert_job(job)
            if (
                job.service_code == "sgnipt"
                and st not in (OrderStatus.COMPLETED, OrderStatus.REPORT_READY)
                and apply_sgnipt_layout_directories(job)
            ):
                self._store.upsert_job(job)
            if st in ACTIVE_BEFORE_RESTART:
                old_log = (job.error_log or "").strip()
                extra = "Daemon restarted while job was active."
                job.error_log = f"{old_log}\n{extra}" if old_log else extra
                job.message = "Interrupted: daemon restarted"
                job.status = OrderStatus.FAILED
                job.completed_at = now_kst_iso()
                self._completed_jobs[job.order_id] = job
                self._stats[job.service_code]["failed"] += 1
                self._store.upsert_job(job)
                interrupted += 1
                logger.warning(
                    "Marked order %s FAILED after restart (was %s)",
                    job.order_id,
                    st.value,
                )
            elif st == OrderStatus.SAVED:
                self._saved_jobs[job.order_id] = job
            elif st == OrderStatus.QUEUED:
                self._jobs[job.order_id] = job
                queued_list.append(job)
                self._stats[job.service_code]["queued"] += 1
            elif st in (OrderStatus.COMPLETED, OrderStatus.REPORT_READY):
                self._completed_jobs[job.order_id] = job
                self._stats[job.service_code]["completed"] += 1
            elif st == OrderStatus.FAILED:
                self._completed_jobs[job.order_id] = job
                self._stats[job.service_code]["failed"] += 1
            elif st == OrderStatus.CANCELLED:
                self._completed_jobs[job.order_id] = job
            else:
                self._completed_jobs[job.order_id] = job
                logger.warning(
                    "Restored order %s with status %s into completed bucket",
                    job.order_id,
                    st.value,
                )

        queued_list.sort(key=lambda j: j.created_at or "")
        for j in queued_list:
            self._queue.put_nowait(j)

        logger.info(
            "Restored %d order(s) from SQLite (queued=%d, interrupted→failed=%d)",
            len(jobs),
            len(queued_list),
            interrupted,
        )

    async def persist_job(self, job: Job) -> None:
        if not self._store:
            return
        try:
            await asyncio.to_thread(self._store.upsert_job, job)
        except Exception as e:
            logger.error("Failed to persist order %s: %s", job.order_id, e)

    async def _ingest_result_snapshot(self, job: Job) -> None:
        if not self._store:
            return
        try:
            await asyncio.to_thread(ingest_result_json_from_disk, self._store, job)
        except Exception as e:
            logger.warning("Ingest result.json for %s: %s", job.order_id, e)

    async def enqueue(self, job: Job) -> int:
        """
        작업을 큐에 추가합니다.

        Args:
            job: 작업 정보

        Returns:
            현재 큐 위치 (1-based)
        """
        if job.service_code == "carrier_screening":
            apply_carrier_layout_directories(job)
        elif job.service_code == "sgnipt":
            apply_sgnipt_layout_directories(job)
        async with self._lock:
            job.status = OrderStatus.QUEUED
            self._jobs[job.order_id] = job
            self._stats[job.service_code]["queued"] += 1

        await self._queue.put(job)
        queue_size = self._queue.qsize()

        logger.info(
            f"[{job.service_code}] Enqueued job {job.order_id} "
            f"(sample: {job.sample_name}, queue_size: {queue_size})"
        )
        await self.persist_job(job)
        return queue_size

    async def dequeue(self) -> Job:
        """큐에서 다음 작업을 가져옵니다 (blocking). 취소 요청된 작업은 건너뜁니다."""
        while True:
            job = await self._queue.get()

            async with self._lock:
                if job.order_id in self._cancel_requested:
                    self._cancel_requested.discard(job.order_id)
                    self._stats[job.service_code]["queued"] = max(
                        0, self._stats[job.service_code]["queued"] - 1
                    )
                    self._jobs.pop(job.order_id, None)
                    job.status = OrderStatus.CANCELLED
                    job.completed_at = now_kst_iso()
                    job.message = "Cancelled while queued"
                    self._completed_jobs[job.order_id] = job
                    cancelled = job
                else:
                    self._stats[job.service_code]["queued"] = max(
                        0, self._stats[job.service_code]["queued"] - 1
                    )
                    cancelled = None

            if cancelled is not None:
                await self.persist_job(cancelled)
                continue

            await self.persist_job(job)
            return job

    async def acquire_slot(self):
        """실행 슬롯 획득 (동시 실행 수 제한)"""
        await self._semaphore.acquire()
        logger.debug(
            f"Acquired execution slot "
            f"(available: {self._semaphore._value}/{self._max_concurrent})"
        )

    def release_slot(self):
        """실행 슬롯 반환"""
        self._semaphore.release()
        logger.debug(
            f"Released execution slot "
            f"(available: {self._semaphore._value}/{self._max_concurrent})"
        )

    async def mark_running(self, job: Job):
        """작업을 실행 중으로 표시"""
        async with self._lock:
            job.status = OrderStatus.RUNNING
            job.started_at = now_kst_iso()
            self._running_jobs[job.order_id] = job
            self._stats[job.service_code]["running"] += 1

        logger.info(f"[{job.service_code}] Job {job.order_id} is now RUNNING")
        await self.persist_job(job)

    async def finalize_reprocess_results(self, job: Job) -> None:
        """
        process_results 재실행 성공 후: result.json 스냅샷 저장, FAILED → COMPLETED 복구.
        COMPLETED / REPORT_READY 는 상태 유지, 메시지·updated_at 만 갱신.
        """
        async with self._lock:
            if job.status == OrderStatus.FAILED:
                self._stats[job.service_code]["failed"] = max(
                    0, self._stats[job.service_code]["failed"] - 1
                )
                self._stats[job.service_code]["completed"] += 1
                job.status = OrderStatus.COMPLETED
                job.completed_at = now_kst_iso()
                job.error_log = None
            job.progress = 100
            job.updated_at = now_kst_iso()
            job.message = "Reprocessed (annotation/QC/result.json)"

        await self._ingest_result_snapshot(job)
        await self.persist_job(job)

    async def mark_completed(self, job: Job):
        """작업을 완료로 표시"""
        async with self._lock:
            job.status = OrderStatus.COMPLETED
            job.completed_at = now_kst_iso()
            self._running_jobs.pop(job.order_id, None)
            self._jobs.pop(job.order_id, None)
            self._completed_jobs[job.order_id] = job
            self._stats[job.service_code]["running"] = max(
                0, self._stats[job.service_code]["running"] - 1
            )
            self._stats[job.service_code]["completed"] += 1

        logger.info(f"[{job.service_code}] Job {job.order_id} COMPLETED")
        await self._ingest_result_snapshot(job)
        await self.persist_job(job)

    async def mark_report_ready(self, order_id: str, message: str = "Report ready for download") -> None:
        """리뷰 후 PDF/HTML 생성 완료 — Portal에서 다운로드 가능."""
        async with self._lock:
            job = (
                self._running_jobs.get(order_id)
                or self._jobs.get(order_id)
                or self._saved_jobs.get(order_id)
                or self._completed_jobs.get(order_id)
            )
            if not job:
                logger.warning("mark_report_ready: order %s not found", order_id)
                return
            job.status = OrderStatus.REPORT_READY
            job.message = message
            job.progress = 100
            job.updated_at = now_kst_iso()
        await self.persist_job(job)
        logger.info("[%s] Order %s → REPORT_READY", job.service_code, order_id)

    async def mark_failed(self, job: Job, error: str = ""):
        """작업을 실패로 표시"""
        async with self._lock:
            job.status = OrderStatus.FAILED
            job.completed_at = now_kst_iso()
            job.error_log = error
            self._running_jobs.pop(job.order_id, None)
            self._jobs.pop(job.order_id, None)
            self._completed_jobs[job.order_id] = job
            self._stats[job.service_code]["running"] = max(
                0, self._stats[job.service_code]["running"] - 1
            )
            self._stats[job.service_code]["failed"] += 1

        logger.error(f"[{job.service_code}] Job {job.order_id} FAILED: {error}")
        await self.persist_job(job)

    def get_job(self, order_id: str) -> Optional[Job]:
        """order_id로 작업 조회"""
        return (
            self._running_jobs.get(order_id)
            or self._jobs.get(order_id)
            or self._saved_jobs.get(order_id)
            or self._completed_jobs.get(order_id)
        )

    async def save_job(self, job: Job) -> None:
        """Portal 저장만: SAVED 상태로 보관, 큐에 넣지 않음."""
        async with self._lock:
            job.status = OrderStatus.SAVED
            job.updated_at = now_kst_iso()
            self._saved_jobs[job.order_id] = job
        if job.service_code == "carrier_screening":
            apply_carrier_layout_directories(job)
        elif job.service_code == "sgnipt":
            apply_sgnipt_layout_directories(job)
        logger.info(f"[{job.service_code}] Saved order {job.order_id} (not queued)")
        await self.persist_job(job)

    async def replace_edited_job(
        self,
        new_job: Job,
        *,
        previous_order_id: Optional[str] = None,
    ) -> None:
        """
        Portal에서 저장·종료·완료 주문 내용을 교체할 때 사용.
        (SAVED 또는 _completed_jobs 의 FAILED / CANCELLED / COMPLETED / REPORT_READY)
        previous_order_id: PATCH URL 등에서 온 기존 키(order_id 변경 시 필수).
        """
        old_id = previous_order_id if previous_order_id is not None else new_job.order_id
        new_id = new_job.order_id

        async with self._lock:
            if old_id not in self._saved_jobs and old_id not in self._completed_jobs:
                raise KeyError(old_id)

            if new_id != old_id:
                for d in (
                    self._saved_jobs,
                    self._jobs,
                    self._running_jobs,
                    self._completed_jobs,
                ):
                    if new_id in d:
                        raise ValueError(f"Order ID already in use: {new_id}")

            if self._store and old_id != new_id:
                try:
                    await asyncio.to_thread(self._store.rename_order_id, old_id, new_id)
                except ValueError:
                    raise
                except Exception as e:
                    logger.error("rename_order_id %s -> %s: %s", old_id, new_id, e)
                    raise ValueError(f"Could not rename order in database: {e}") from e

            if old_id in self._saved_jobs:
                if self._saved_jobs[old_id].status != OrderStatus.SAVED:
                    raise ValueError(f"Order {old_id} is not SAVED")
                self._saved_jobs.pop(old_id)
                self._saved_jobs[new_id] = new_job
            else:
                st = self._completed_jobs[old_id].status
                if st not in (
                    OrderStatus.FAILED,
                    OrderStatus.CANCELLED,
                    OrderStatus.COMPLETED,
                    OrderStatus.REPORT_READY,
                ):
                    raise ValueError(
                        f"Order {old_id} cannot be edited in status {st.value}"
                    )
                self._completed_jobs.pop(old_id)
                self._completed_jobs[new_id] = new_job

        if new_job.service_code == "carrier_screening":
            apply_carrier_layout_directories(new_job)
        elif new_job.service_code == "sgnipt":
            apply_sgnipt_layout_directories(new_job)
        logger.info(
            "[%s] Updated order %s -> %s (status=%s)",
            new_job.service_code,
            old_id,
            new_id,
            new_job.status.value,
        )
        await self.persist_job(new_job)

    def _prepare_job_for_retry(self, job: Job) -> None:
        """FAILED/CANCELLED/COMPLETED 주문을 다시 큐에 넣기 전 필드·통계 정리."""
        if job.status == OrderStatus.FAILED:
            self._stats[job.service_code]["failed"] = max(
                0, self._stats[job.service_code]["failed"] - 1
            )
        elif job.status in (OrderStatus.COMPLETED, OrderStatus.REPORT_READY):
            self._stats[job.service_code]["completed"] = max(
                0, self._stats[job.service_code]["completed"] - 1
            )
        job.progress = 0
        job.message = ""
        job.error_log = None
        job.started_at = None
        job.completed_at = None
        job.updated_at = now_kst_iso()
        job.pid = None
        job.exit_code = None
        job.duration = None

    async def start_saved_job(self, order_id: str) -> tuple:
        """
        SAVED 주문, 또는 FAILED/CANCELLED/COMPLETED 주문(재분석)을 큐에 넣습니다.
        (Job, queue_position) 반환.
        """
        async with self._lock:
            job = self._saved_jobs.pop(order_id, None)
            if job:
                if job.status != OrderStatus.SAVED:
                    self._saved_jobs[order_id] = job
                    raise ValueError(f"Order {order_id} is not SAVED")
            else:
                cj = self._completed_jobs.get(order_id)
                if cj is None or cj.status not in (
                    OrderStatus.FAILED,
                    OrderStatus.CANCELLED,
                    OrderStatus.COMPLETED,
                    OrderStatus.REPORT_READY,
                ):
                    raise KeyError(order_id)
                job = self._completed_jobs.pop(order_id)
                self._prepare_job_for_retry(job)

        queue_position = await self.enqueue(job)
        return job, queue_position

    async def delete_saved(self, order_id: str) -> bool:
        """SAVED 주문만 삭제."""
        async with self._lock:
            if order_id not in self._saved_jobs:
                return False
            del self._saved_jobs[order_id]
        if self._store:
            try:
                await asyncio.to_thread(self._store.delete_order, order_id)
            except Exception as e:
                logger.error("Failed to delete order %s from store: %s", order_id, e)
        logger.info(f"Deleted saved order {order_id}")
        return True

    async def forget_order(self, order_id: str) -> Optional[Job]:
        """
        order_id를 저장/큐/실행/완료 버킷과 asyncio.Queue에서 제거합니다 (SQLite 제외).
        통계 카운터를 Job 상태에 맞게 감소시킵니다.
        """
        async with self._lock:
            job = (
                self._saved_jobs.pop(order_id, None)
                or self._jobs.pop(order_id, None)
                or self._running_jobs.pop(order_id, None)
                or self._completed_jobs.pop(order_id, None)
            )
            self._saved_jobs.pop(order_id, None)
            self._jobs.pop(order_id, None)
            self._running_jobs.pop(order_id, None)
            self._completed_jobs.pop(order_id, None)
            self._cancel_requested.discard(order_id)
            pending: List[Job] = []
            while True:
                try:
                    j = self._queue.get_nowait()
                except asyncio.QueueEmpty:
                    break
                if j.order_id != order_id:
                    pending.append(j)
            for j in pending:
                self._queue.put_nowait(j)
            if job:
                svc = job.service_code
                st = job.status
                if st == OrderStatus.QUEUED:
                    self._stats[svc]["queued"] = max(
                        0, self._stats[svc]["queued"] - 1
                    )
                elif st in (
                    OrderStatus.RUNNING,
                    OrderStatus.DOWNLOADING,
                    OrderStatus.PROCESSING,
                    OrderStatus.UPLOADING,
                    OrderStatus.RECEIVED,
                ):
                    self._stats[svc]["running"] = max(
                        0, self._stats[svc]["running"] - 1
                    )
                elif st in (OrderStatus.COMPLETED, OrderStatus.REPORT_READY):
                    self._stats[svc]["completed"] = max(
                        0, self._stats[svc]["completed"] - 1
                    )
                elif st == OrderStatus.FAILED:
                    self._stats[svc]["failed"] = max(
                        0, self._stats[svc]["failed"] - 1
                    )
        return job

    async def purge_queued_order(self, order_id: str) -> Optional[Job]:
        """
        QUEUED 주문을 메모리 큐와 asyncio.Queue에서 제거합니다.
        반환: 제거된 Job (없으면 None).
        """
        async with self._lock:
            job = self._jobs.pop(order_id, None)
            if not job:
                return None
            self._stats[job.service_code]["queued"] = max(
                0, self._stats[job.service_code]["queued"] - 1
            )
            self._cancel_requested.discard(order_id)
            pending: List[Job] = []
            while True:
                try:
                    j = self._queue.get_nowait()
                except asyncio.QueueEmpty:
                    break
                if j.order_id != order_id:
                    pending.append(j)
            for j in pending:
                self._queue.put_nowait(j)
        logger.info("Purged queued order %s from asyncio queue", order_id)
        return job

    async def delete_order_with_artifacts(self, order_id: str) -> tuple:
        """
        주문 레코드를 제거하고 analysis/output/log 트리 삭제 (FASTQ 제외).
        Returns:
            (ok: bool, message: str, detail: dict)
        """
        job = self.get_job(order_id)
        if not job:
            return False, f"Order not found: {order_id}", {}

        blocked = (
            OrderStatus.RUNNING,
            OrderStatus.DOWNLOADING,
            OrderStatus.PROCESSING,
            OrderStatus.UPLOADING,
            OrderStatus.RECEIVED,
        )
        if job.status in blocked:
            return (
                False,
                "Cannot delete while the job is active; use Stop first, then Delete.",
                {},
            )

        if job.service_code == "carrier_screening":
            apply_carrier_layout_directories(job)
        elif job.service_code == "sgnipt":
            apply_sgnipt_layout_directories(job)

        if job.status == OrderStatus.QUEUED:
            pj = await self.purge_queued_order(order_id)
            if pj is None:
                logger.warning(
                    "delete_order_with_artifacts: QUEUED %s not in asyncio queue; forcing forget",
                    order_id,
                )
                await self.forget_order(order_id)
            else:
                job = pj

        deleted, errs = await asyncio.to_thread(delete_run_artifacts, job)

        forgotten = await self.forget_order(order_id)

        if self._store:
            try:
                await asyncio.to_thread(self._store.delete_order, order_id)
            except Exception as e:
                logger.error("delete_order from store for %s: %s", order_id, e)

        detail = {"deleted": deleted, "errors": errs, "memory_cleared": forgotten is not None}
        msg = f"Order {order_id} removed; deleted {len(deleted)} directory tree(s)"
        if errs:
            msg += f"; {len(errs)} path warning(s)"
        return True, msg, detail

    async def purge_order_db_only(self, order_id: str, force: bool = False) -> tuple:
        """
        SQLite 행과 메모리·큐에서만 제거합니다 (디스크 analysis/output/log 미삭제).
        job_json 손상으로 부팅 시 복원되지 않은 행도 SQLite에서 지울 수 있습니다.
        force=True 이면 RUNNING 등 활성 상태여도 DB/메모리에서 제거합니다(파이프라인 프로세스는 그대로).
        """
        job = self.get_job(order_id)
        blocked = (
            OrderStatus.RUNNING,
            OrderStatus.DOWNLOADING,
            OrderStatus.PROCESSING,
            OrderStatus.UPLOADING,
            OrderStatus.RECEIVED,
        )
        if job and job.status in blocked and not force:
            return (
                False,
                "Stop the job before removing this record from the database, or use Purge with Force.",
                {},
            )
        if job and job.status in blocked and force:
            logger.warning(
                "purge_order_db_only force=True for active job %s status=%s",
                order_id,
                job.status.value,
            )

        forgotten = await self.forget_order(order_id)
        sqlite_rows = 0
        if self._store:
            try:
                sqlite_rows = await asyncio.to_thread(
                    self._store.delete_order, order_id
                )
            except Exception as e:
                logger.error("purge_order_db_only SQLite %s: %s", order_id, e)
                return False, f"SQLite delete failed: {e}", {}

        msg = f"Order {order_id} removed from database and daemon memory"
        if forgotten is None and job is None:
            msg += " (no live job in memory; SQLite row deleted if it existed)"
        detail = {
            "had_memory_job": forgotten is not None,
            "sqlite_rows_deleted": sqlite_rows,
        }
        if sqlite_rows == 0:
            msg += " — SQLite: no row matched this order_id (already removed or typo)."
        return True, msg, detail

    async def request_cancel_queued(self, order_id: str) -> bool:
        """QUEUED 상태 주문에 대해 dequeue 시 버리도록 표시."""
        async with self._lock:
            job = self._jobs.get(order_id)
            if not job or job.status != OrderStatus.QUEUED:
                return False
            self._cancel_requested.add(order_id)
        logger.info(f"Cancel requested for queued order {order_id}")
        return True

    async def mark_cancelled(self, job: Job, message: str = "Cancelled by user"):
        """실행 중 사용자 취소 등 (실패 통계 없이 종료 처리)."""
        async with self._lock:
            job.status = OrderStatus.CANCELLED
            job.completed_at = now_kst_iso()
            job.message = message
            self._running_jobs.pop(job.order_id, None)
            self._jobs.pop(job.order_id, None)
            self._saved_jobs.pop(job.order_id, None)
            self._completed_jobs[job.order_id] = job
            self._stats[job.service_code]["running"] = max(
                0, self._stats[job.service_code]["running"] - 1
            )
        logger.info(f"[{job.service_code}] Job {job.order_id} CANCELLED")
        await self.persist_job(job)

    async def update_queued_job_fastq_paths(
        self,
        order_id: str,
        fastq_r1_path: Optional[str] = None,
        fastq_r2_path: Optional[str] = None,
    ) -> tuple[bool, str]:
        """
        QUEUED 상태인 주문만 로컬 FASTQ 경로 변경.
        fastq_r*_path 가 None 이면 해당 필드는 변경하지 않음.
        """
        async with self._lock:
            job = self._jobs.get(order_id) or self._saved_jobs.get(order_id)
            if not job:
                return False, "Order not found"
            if job.status not in (OrderStatus.QUEUED, OrderStatus.SAVED):
                return False, "Only SAVED or QUEUED orders can update FASTQ paths"
            if fastq_r1_path is not None:
                job.fastq_r1_path = fastq_r1_path or None
            if fastq_r2_path is not None:
                job.fastq_r2_path = fastq_r2_path or None
            job.updated_at = now_kst_iso()
        if job.service_code == "carrier_screening":
            apply_carrier_layout_directories(job)
        elif job.service_code == "sgnipt":
            apply_sgnipt_layout_directories(job)
        await self.persist_job(job)
        return True, ""

    def get_running_jobs(self, service_code: Optional[str] = None) -> List[Job]:
        """실행 중인 작업 목록"""
        jobs = list(self._running_jobs.values())
        if service_code:
            jobs = [j for j in jobs if j.service_code == service_code]
        return jobs

    def get_summary(self) -> QueueSummary:
        """큐 상태 요약"""
        total_queued = self._queue.qsize()
        total_running = len(self._running_jobs)
        total_completed = sum(s["completed"] for s in self._stats.values())
        total_failed = sum(s["failed"] for s in self._stats.values())

        jobs_by_service = {}
        for svc, stats in self._stats.items():
            jobs_by_service[svc] = stats["queued"] + stats["running"]

        running_jobs = [
            {
                "order_id": j.order_id,
                "service_code": j.service_code,
                "sample_name": j.sample_name,
                "status": j.status.value,
                "progress": j.progress,
                "started_at": j.started_at,
            }
            for j in self._running_jobs.values()
        ]

        return QueueSummary(
            total_queued=total_queued,
            total_running=total_running,
            total_completed=total_completed,
            total_failed=total_failed,
            jobs_by_service=jobs_by_service,
            running_jobs=running_jobs,
        )

    @property
    def available_slots(self) -> int:
        """사용 가능한 실행 슬롯 수"""
        return self._semaphore._value

    @property
    def max_concurrent(self) -> int:
        return self._max_concurrent


# 전역 인스턴스
_queue_manager: Optional[QueueManager] = None


def get_queue_manager() -> QueueManager:
    """전역 QueueManager 인스턴스 반환"""
    global _queue_manager
    if _queue_manager is None:
        store = OrderStore(settings.resolved_orders_db_path)
        _queue_manager = QueueManager(store=store)
    return _queue_manager
