"""
Queue Manager

작업 큐와 동시 실행 제어를 관리합니다.
기존 nipt-daemon의 queue_manager.py를 일반화하여 멀티서비스를 지원합니다.
"""

import asyncio
import logging
from typing import Dict, List, Optional, Any
from datetime import datetime
from collections import defaultdict

from .config import settings
from .models import Job, OrderStatus, QueueSummary

logger = logging.getLogger(__name__)


class QueueManager:
    """
    멀티서비스 작업 큐 관리자.
    
    단일 큐에서 모든 서비스의 작업을 관리하며,
    Semaphore를 통해 동시 실행 수를 제어합니다.
    """

    def __init__(self, max_concurrent: int = None):
        self._max_concurrent = max_concurrent or settings.max_concurrent_jobs
        self._queue: asyncio.Queue = asyncio.Queue()
        self._semaphore = asyncio.Semaphore(self._max_concurrent)
        
        # 작업 상태 추적
        self._jobs: Dict[str, Job] = {}           # order_id -> Job
        self._running_jobs: Dict[str, Job] = {}    # order_id -> Job (실행 중)
        self._completed_jobs: Dict[str, Job] = {}  # order_id -> Job (완료/실패)
        
        # 서비스별 통계
        self._stats: Dict[str, Dict[str, int]] = defaultdict(
            lambda: {"queued": 0, "running": 0, "completed": 0, "failed": 0}
        )
        
        self._lock = asyncio.Lock()
        logger.info(f"QueueManager initialized (max_concurrent={self._max_concurrent})")

    async def enqueue(self, job: Job) -> int:
        """
        작업을 큐에 추가합니다.
        
        Args:
            job: 작업 정보
            
        Returns:
            현재 큐 위치 (1-based)
        """
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
        return queue_size

    async def dequeue(self) -> Job:
        """큐에서 다음 작업을 가져옵니다 (blocking)."""
        job = await self._queue.get()
        
        async with self._lock:
            self._stats[job.service_code]["queued"] = max(
                0, self._stats[job.service_code]["queued"] - 1
            )
        
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
            job.started_at = datetime.now().isoformat()
            self._running_jobs[job.order_id] = job
            self._stats[job.service_code]["running"] += 1
        
        logger.info(f"[{job.service_code}] Job {job.order_id} is now RUNNING")

    async def mark_completed(self, job: Job):
        """작업을 완료로 표시"""
        async with self._lock:
            job.status = OrderStatus.COMPLETED
            job.completed_at = datetime.now().isoformat()
            self._running_jobs.pop(job.order_id, None)
            self._completed_jobs[job.order_id] = job
            self._stats[job.service_code]["running"] = max(
                0, self._stats[job.service_code]["running"] - 1
            )
            self._stats[job.service_code]["completed"] += 1
        
        logger.info(f"[{job.service_code}] Job {job.order_id} COMPLETED")

    async def mark_failed(self, job: Job, error: str = ""):
        """작업을 실패로 표시"""
        async with self._lock:
            job.status = OrderStatus.FAILED
            job.completed_at = datetime.now().isoformat()
            job.error_log = error
            self._running_jobs.pop(job.order_id, None)
            self._completed_jobs[job.order_id] = job
            self._stats[job.service_code]["running"] = max(
                0, self._stats[job.service_code]["running"] - 1
            )
            self._stats[job.service_code]["failed"] += 1
        
        logger.error(f"[{job.service_code}] Job {job.order_id} FAILED: {error}")

    def get_job(self, order_id: str) -> Optional[Job]:
        """order_id로 작업 조회"""
        return (
            self._running_jobs.get(order_id) or
            self._jobs.get(order_id) or
            self._completed_jobs.get(order_id)
        )

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
        total_completed = sum(
            s["completed"] for s in self._stats.values()
        )
        total_failed = sum(
            s["failed"] for s in self._stats.values()
        )
        
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
                "started_at": j.started_at
            }
            for j in self._running_jobs.values()
        ]

        return QueueSummary(
            total_queued=total_queued,
            total_running=total_running,
            total_completed=total_completed,
            total_failed=total_failed,
            jobs_by_service=jobs_by_service,
            running_jobs=running_jobs
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
        _queue_manager = QueueManager()
    return _queue_manager
