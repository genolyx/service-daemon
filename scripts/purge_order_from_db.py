#!/usr/bin/env python3
"""SQLite orders.db 에서 주문 행만 삭제 (데몬 없이)."""

import argparse
import os
import sqlite3
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


def main() -> None:
    parser = argparse.ArgumentParser(description="DELETE FROM orders WHERE order_id = …")
    parser.add_argument("order_id", help="order_id PK (예: NA12878_test_rebuild)")
    parser.add_argument(
        "--db",
        default=None,
        help="orders.db 경로 (기본: settings.resolved_orders_db_path)",
    )
    args = parser.parse_args()

    if args.db:
        db_path = os.path.abspath(args.db)
    else:
        from app.config import settings

        db_path = settings.resolved_orders_db_path

    if not os.path.isfile(db_path):
        print(f"DB not found: {db_path}", file=sys.stderr)
        sys.exit(1)

    conn = sqlite3.connect(db_path)
    try:
        cur = conn.execute("DELETE FROM orders WHERE order_id = ?", (args.order_id,))
        conn.commit()
        print(f"Deleted {cur.rowcount} row(s) for order_id={args.order_id!r} in {db_path}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
