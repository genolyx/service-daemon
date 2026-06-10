"""
Portal variant sets: user-uploaded TSV lists tagged by name (e.g. Hotspot).

TSV columns: chrom, pos, ref, alt (required); gene, label (optional).
Comment lines (#) and blank lines are skipped.
"""

from __future__ import annotations

import os
import sqlite3
import threading
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from .config import settings

_LOCK = threading.Lock()
_REQUIRED_COLS = frozenset({"chrom", "pos", "ref", "alt"})


def resolved_variant_sets_db_path() -> str:
    return os.path.join(settings.base_dir, "portal", "variant_sets.db")


def norm_chrom(chrom: str) -> str:
    c = (chrom or "").strip().lower().replace(" ", "")
    if not c:
        return ""
    return c if c.startswith("chr") else f"chr{c}"


def variant_lookup_key(chrom: str, pos: Any, ref: str, alt: str) -> Optional[str]:
    nc = norm_chrom(chrom)
    if not nc:
        return None
    try:
        ipos = int(str(pos).strip())
    except (TypeError, ValueError):
        return None
    r = (ref or "").strip().upper()
    a = (alt or "").strip().upper()
    if not r or not a:
        return None
    return f"{nc}:{ipos}:{r}:{a}"


def _connect() -> sqlite3.Connection:
    path = resolved_variant_sets_db_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    conn = sqlite3.connect(path, timeout=30)
    conn.row_factory = sqlite3.Row
    return conn


def init_variant_sets_db() -> None:
    with _LOCK:
        conn = _connect()
        try:
            conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS variant_sets (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    tag_name TEXT NOT NULL UNIQUE COLLATE NOCASE,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS variant_set_entries (
                    set_id INTEGER NOT NULL,
                    chrom TEXT NOT NULL,
                    pos INTEGER NOT NULL,
                    ref TEXT NOT NULL,
                    alt TEXT NOT NULL,
                    gene TEXT,
                    label TEXT,
                    FOREIGN KEY (set_id) REFERENCES variant_sets(id) ON DELETE CASCADE
                );
                CREATE INDEX IF NOT EXISTS idx_vse_set ON variant_set_entries(set_id);
                CREATE INDEX IF NOT EXISTS idx_vse_locus ON variant_set_entries(chrom, pos, ref, alt);
                """
            )
            conn.commit()
        finally:
            conn.close()


def _now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def parse_variant_sets_tsv(text: str) -> List[Dict[str, Any]]:
    """Parse TSV text into entry dicts. Raises ValueError on invalid format."""
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    header: Optional[List[str]] = None
    entries: List[Dict[str, Any]] = []
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            maybe = line.lstrip("#").strip()
            maybe_parts = maybe.split("\t")
            maybe_hdr = [h.strip().lower() for h in maybe_parts if h.strip()]
            if maybe_hdr and _REQUIRED_COLS.issubset(set(maybe_hdr)):
                header = maybe_hdr
            continue
        parts = line.split("\t")
        if header is None:
            header = [h.strip().lower() for h in parts]
            missing = _REQUIRED_COLS - set(header)
            if missing:
                raise ValueError(
                    f"TSV header must include: chrom, pos, ref, alt (missing: {', '.join(sorted(missing))})"
                )
            continue
        if len(parts) < len(header):
            parts.extend([""] * (len(header) - len(parts)))
        row = {header[i]: (parts[i].strip() if i < len(parts) else "") for i in range(len(header))}
        key = variant_lookup_key(row.get("chrom", ""), row.get("pos", ""), row.get("ref", ""), row.get("alt", ""))
        if not key:
            continue
        entries.append(
            {
                "chrom": norm_chrom(row.get("chrom", "")),
                "pos": int(str(row.get("pos", "")).strip()),
                "ref": row.get("ref", "").strip().upper(),
                "alt": row.get("alt", "").strip().upper(),
                "gene": (row.get("gene") or "").strip() or None,
                "label": (row.get("label") or "").strip() or None,
            }
        )
    if header is None:
        raise ValueError("TSV has no header row (expected: chrom\\tpos\\tref\\talt...)")
    if not entries:
        raise ValueError("TSV contains no valid variant rows")
    return entries


def list_variant_sets() -> List[Dict[str, Any]]:
    init_variant_sets_db()
    with _LOCK:
        conn = _connect()
        try:
            rows = conn.execute(
                """
                SELECT s.id, s.tag_name, s.created_at, s.updated_at,
                       COUNT(e.rowid) AS entry_count
                FROM variant_sets s
                LEFT JOIN variant_set_entries e ON e.set_id = s.id
                GROUP BY s.id
                ORDER BY s.tag_name COLLATE NOCASE
                """
            ).fetchall()
            return [dict(r) for r in rows]
        finally:
            conn.close()


def export_lookup() -> Dict[str, List[str]]:
    """Map 'chr:pos:REF:ALT' -> sorted tag names."""
    init_variant_sets_db()
    lookup: Dict[str, List[str]] = {}
    with _LOCK:
        conn = _connect()
        try:
            rows = conn.execute(
                """
                SELECT s.tag_name, e.chrom, e.pos, e.ref, e.alt
                FROM variant_set_entries e
                JOIN variant_sets s ON s.id = e.set_id
                """
            ).fetchall()
            for r in rows:
                key = variant_lookup_key(r["chrom"], r["pos"], r["ref"], r["alt"])
                if not key:
                    continue
                lookup.setdefault(key, [])
                tag = r["tag_name"]
                if tag not in lookup[key]:
                    lookup[key].append(tag)
            for key in lookup:
                lookup[key].sort(key=str.lower)
            return lookup
        finally:
            conn.close()


def export_entries_by_tag() -> Dict[str, List[Dict[str, Any]]]:
    """All variants grouped by tag name (for Portal expand panel)."""
    init_variant_sets_db()
    by_tag: Dict[str, List[Dict[str, Any]]] = {}
    with _LOCK:
        conn = _connect()
        try:
            rows = conn.execute(
                """
                SELECT s.tag_name, e.chrom, e.pos, e.ref, e.alt, e.gene, e.label
                FROM variant_set_entries e
                JOIN variant_sets s ON s.id = e.set_id
                ORDER BY s.tag_name COLLATE NOCASE, e.chrom, e.pos, e.ref, e.alt
                """
            ).fetchall()
            for r in rows:
                tag = r["tag_name"]
                by_tag.setdefault(tag, []).append(
                    {
                        "chrom": r["chrom"],
                        "pos": int(r["pos"]),
                        "ref": r["ref"],
                        "alt": r["alt"],
                        "gene": r["gene"],
                        "label": r["label"],
                    }
                )
            return by_tag
        finally:
            conn.close()


def get_catalog_for_portal() -> Dict[str, Any]:
    return {
        "sets": list_variant_sets(),
        "lookup": export_lookup(),
        "entries_by_tag": export_entries_by_tag(),
    }


def upsert_variant_set(tag_name: str, entries: List[Dict[str, Any]]) -> Dict[str, Any]:
    tag = (tag_name or "").strip()
    if not tag:
        raise ValueError("tag_name is required")
    if len(tag) > 80:
        raise ValueError("tag_name must be 80 characters or fewer")
    init_variant_sets_db()
    now = _now_iso()
    with _LOCK:
        conn = _connect()
        try:
            conn.execute("PRAGMA foreign_keys = ON")
            existing = conn.execute(
                "SELECT id FROM variant_sets WHERE tag_name = ? COLLATE NOCASE", (tag,)
            ).fetchone()
            if existing:
                set_id = int(existing["id"])
                conn.execute("DELETE FROM variant_set_entries WHERE set_id = ?", (set_id,))
                conn.execute(
                    "UPDATE variant_sets SET updated_at = ? WHERE id = ?", (now, set_id)
                )
            else:
                cur = conn.execute(
                    "INSERT INTO variant_sets (tag_name, created_at, updated_at) VALUES (?, ?, ?)",
                    (tag, now, now),
                )
                set_id = int(cur.lastrowid)
            conn.executemany(
                """
                INSERT INTO variant_set_entries (set_id, chrom, pos, ref, alt, gene, label)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        set_id,
                        e["chrom"],
                        e["pos"],
                        e["ref"],
                        e["alt"],
                        e.get("gene"),
                        e.get("label"),
                    )
                    for e in entries
                ],
            )
            conn.commit()
            return {
                "id": set_id,
                "tag_name": tag,
                "entry_count": len(entries),
                "updated_at": now,
            }
        finally:
            conn.close()


def delete_variant_set(set_id: int) -> bool:
    init_variant_sets_db()
    with _LOCK:
        conn = _connect()
        try:
            conn.execute("PRAGMA foreign_keys = ON")
            cur = conn.execute("DELETE FROM variant_sets WHERE id = ?", (set_id,))
            conn.commit()
            return cur.rowcount > 0
        finally:
            conn.close()
