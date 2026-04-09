"""
SQLite persistence for orders (Job snapshots + result/review/report JSON blobs).

Default path: {BASE_DIR}/service-daemon/orders.db — mount BASE_DIR in Docker for durability.
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
import threading
from typing import Any, Dict, List, Optional

from .datetime_kst import now_kst_iso
from .models import Job, OrderStatus

logger = logging.getLogger(__name__)


class OrderStore:
    """Thread-safe SQLite store for Job rows and optional JSON documents per order."""

    def __init__(self, db_path: str):
        self._path = os.path.abspath(db_path)
        parent = os.path.dirname(self._path)
        if parent:
            os.makedirs(parent, mode=0o755, exist_ok=True)
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(self._path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        with self._lock:
            self._conn.execute("PRAGMA journal_mode=WAL")
        self._init_schema()

    def _init_schema(self) -> None:
        with self._lock:
            self._conn.execute(
                """
                CREATE TABLE IF NOT EXISTS orders (
                    order_id TEXT PRIMARY KEY NOT NULL,
                    status TEXT NOT NULL,
                    job_json TEXT NOT NULL,
                    result_json TEXT,
                    review_json TEXT,
                    report_json TEXT,
                    extra_json TEXT,
                    updated_at TEXT NOT NULL
                )
                """
            )
            self._conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_orders_status ON orders(status)"
            )
            self._conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_orders_updated ON orders(updated_at)"
            )
            self._conn.commit()

    def close(self) -> None:
        with self._lock:
            self._conn.close()

    def delete_order(self, order_id: str) -> int:
        with self._lock:
            cur = self._conn.execute(
                "DELETE FROM orders WHERE order_id = ?", (order_id,)
            )
            self._conn.commit()
            return int(cur.rowcount or 0)

    def rename_order_id(self, old_id: str, new_id: str) -> None:
        """PK 이전: result/review/report 등 기존 열 유지. old 행이 없으면 무시."""
        if old_id == new_id:
            return
        now = now_kst_iso()
        with self._lock:
            cur = self._conn.cursor()
            cur.execute(
                """
                SELECT status, job_json, result_json, review_json, report_json, extra_json
                FROM orders WHERE order_id = ?
                """,
                (old_id,),
            )
            row = cur.fetchone()
            if not row:
                return
            cur.execute("SELECT 1 FROM orders WHERE order_id = ?", (new_id,))
            if cur.fetchone():
                raise ValueError(f"Order ID already in database: {new_id}")
            cur.execute(
                """
                INSERT INTO orders (
                    order_id, status, job_json, result_json, review_json,
                    report_json, extra_json, updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (new_id, row[0], row[1], row[2], row[3], row[4], row[5], now),
            )
            cur.execute("DELETE FROM orders WHERE order_id = ?", (old_id,))
            self._conn.commit()

    def upsert_job(self, job: Job) -> None:
        """Write full Job JSON; preserve existing document columns."""
        oid = job.order_id
        job_json = json.dumps(job.model_dump(mode="json"), ensure_ascii=False)
        now = now_kst_iso()
        with self._lock:
            cur = self._conn.cursor()
            cur.execute(
                "SELECT result_json, review_json, report_json, extra_json FROM orders WHERE order_id = ?",
                (oid,),
            )
            row = cur.fetchone()
            rj, rvw, rpt, exj = (None, None, None, None)
            if row:
                rj, rvw, rpt, exj = row[0], row[1], row[2], row[3]
            cur.execute(
                """
                INSERT INTO orders (order_id, status, job_json, result_json, review_json, report_json, extra_json, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(order_id) DO UPDATE SET
                    status = excluded.status,
                    job_json = excluded.job_json,
                    updated_at = excluded.updated_at
                """,
                (oid, job.status.value, job_json, rj, rvw, rpt, exj, now),
            )
            self._conn.commit()

    def set_result_json(self, order_id: str, data: Dict[str, Any]) -> None:
        s = json.dumps(data, ensure_ascii=False)
        now = now_kst_iso()
        with self._lock:
            self._conn.execute(
                "UPDATE orders SET result_json = ?, updated_at = ? WHERE order_id = ?",
                (s, now, order_id),
            )
            self._conn.commit()

    def clear_result_json(self, order_id: str) -> None:
        """Remove cached result_json (e.g. stale snapshot after fixing per-order disk paths)."""
        now = now_kst_iso()
        with self._lock:
            self._conn.execute(
                "UPDATE orders SET result_json = NULL, updated_at = ? WHERE order_id = ?",
                (now, order_id),
            )
            self._conn.commit()

    def set_review_json(self, order_id: str, data: Dict[str, Any]) -> None:
        s = json.dumps(data, ensure_ascii=False)
        now = now_kst_iso()
        with self._lock:
            self._conn.execute(
                "UPDATE orders SET review_json = ?, updated_at = ? WHERE order_id = ?",
                (s, now, order_id),
            )
            self._conn.commit()

    def set_report_json(self, order_id: str, data: Dict[str, Any]) -> None:
        s = json.dumps(data, ensure_ascii=False)
        now = now_kst_iso()
        with self._lock:
            self._conn.execute(
                "UPDATE orders SET report_json = ?, updated_at = ? WHERE order_id = ?",
                (s, now, order_id),
            )
            self._conn.commit()

    def set_extra_json(self, order_id: str, data: Dict[str, Any]) -> None:
        s = json.dumps(data, ensure_ascii=False)
        now = now_kst_iso()
        with self._lock:
            self._conn.execute(
                "UPDATE orders SET extra_json = ?, updated_at = ? WHERE order_id = ?",
                (s, now, order_id),
            )
            self._conn.commit()

    def get_result_json(self, order_id: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            cur = self._conn.cursor()
            cur.execute(
                "SELECT result_json FROM orders WHERE order_id = ?", (order_id,)
            )
            row = cur.fetchone()
            if not row or not row[0]:
                return None
            try:
                return json.loads(row[0])
            except json.JSONDecodeError:
                return None

    def get_review_json(self, order_id: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            cur = self._conn.cursor()
            cur.execute(
                "SELECT review_json FROM orders WHERE order_id = ?", (order_id,)
            )
            row = cur.fetchone()
            if not row or not row[0]:
                return None
            try:
                return json.loads(row[0])
            except json.JSONDecodeError:
                return None

    def get_report_json(self, order_id: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            cur = self._conn.cursor()
            cur.execute(
                "SELECT report_json FROM orders WHERE order_id = ?", (order_id,)
            )
            row = cur.fetchone()
            if not row or not row[0]:
                return None
            try:
                return json.loads(row[0])
            except json.JSONDecodeError:
                return None

    def fetch_job(self, order_id: str) -> Optional[Job]:
        """Load a single Job snapshot by primary key (for prior-order resolution, etc.)."""
        with self._lock:
            cur = self._conn.cursor()
            cur.execute("SELECT job_json FROM orders WHERE order_id = ?", (order_id,))
            row = cur.fetchone()
        if not row or not row[0]:
            return None
        try:
            return Job.model_validate(json.loads(row[0]))
        except Exception as e:
            logger.warning("fetch_job %s: %s", order_id, e)
            return None

    def fetch_all_jobs(self, *, exclude_statuses: Optional[List[str]] = None) -> List[Job]:
        """Load Job rows from job_json, optionally excluding certain statuses."""
        with self._lock:
            cur = self._conn.cursor()
            if exclude_statuses:
                placeholders = ",".join("?" for _ in exclude_statuses)
                cur.execute(
                    f"SELECT job_json FROM orders WHERE status NOT IN ({placeholders})",
                    exclude_statuses,
                )
            else:
                cur.execute("SELECT job_json FROM orders")
            rows = cur.fetchall()
        out: List[Job] = []
        for r in rows:
            if not r[0]:
                continue
            try:
                out.append(Job.model_validate(json.loads(r[0])))
            except Exception as e:
                logger.warning("Skipping corrupt job row: %s", e)
        return out


def ingest_result_json_from_disk(store: OrderStore, job: Job) -> None:
    """If result.json exists for this order (carrier: canonical review path), store in DB."""
    path = None
    _CARRIER_LIKE = {"carrier_screening", "whole_exome", "health_screening"}
    if getattr(job, "service_code", None) in _CARRIER_LIKE:
        from app.services.carrier_screening.plugin import carrier_result_json_path

        path = carrier_result_json_path(job)
    elif job.output_dir:
        path = os.path.join(job.output_dir, "result.json")
    if not path or not os.path.isfile(path):
        return
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            store.set_result_json(job.order_id, data)
    except Exception as e:
        logger.warning("Could not ingest result.json for %s: %s", job.order_id, e)


def ingest_report_json_from_disk(
    store: OrderStore, job: Job, output_dir: Optional[str] = None
) -> None:
    """If output_dir/report.json exists, store in DB."""
    od = (output_dir or "").strip() or job.output_dir
    if not od:
        return
    path = os.path.join(od, "report.json")
    if not os.path.isfile(path):
        return
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            store.set_report_json(job.order_id, data)
    except Exception as e:
        logger.warning("Could not ingest report.json for %s: %s", job.order_id, e)


ACTIVE_BEFORE_RESTART = frozenset(
    {
        OrderStatus.RUNNING,
        OrderStatus.DOWNLOADING,
        OrderStatus.PROCESSING,
        OrderStatus.UPLOADING,
        OrderStatus.RECEIVED,
    }
)
