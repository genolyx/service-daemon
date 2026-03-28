#!/usr/bin/env python3
"""
Batch-fill ``gene_data`` SQLite (Sam-compatible schema) using Gemini + Google Search.

Usage:
  export GEMINI_API_KEY=...
  python scripts/populate_gene_knowledge_db.py --db /path/to/gene_knowledge.db CFTR SMN1 ALPL

Or genes from stdin (one per line):
  echo -e "CFTR\\nSMN1" | python scripts/populate_gene_knowledge_db.py --db ./data/db/gene_knowledge.db -

Requires: google-genai, pydantic (see requirements.txt).
"""
from __future__ import annotations

import argparse
import os
import sys

# Repo root on path
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)


def main() -> int:
    p = argparse.ArgumentParser(description="Populate gene_knowledge SQLite via Gemini")
    p.add_argument(
        "--db",
        required=True,
        help="Path to SQLite file (created if missing)",
    )
    p.add_argument(
        "--model",
        default=os.environ.get("GENE_KNOWLEDGE_GEMINI_MODEL", "gemini-2.5-flash"),
    )
    p.add_argument(
        "genes",
        nargs="*",
        help="Gene symbols; use - to read from stdin",
    )
    args = p.parse_args()

    api_key = os.environ.get("GEMINI_API_KEY") or os.environ.get("GEMINI_API_KEY".lower())
    if not api_key:
        print("Set GEMINI_API_KEY", file=sys.stderr)
        return 1

    from app.services.carrier_screening.gene_knowledge_db import (
        get_or_fetch_gene_diseases,
        init_gene_knowledge_database,
    )

    init_gene_knowledge_database(args.db)

    genes: list[str] = []
    if args.genes == ["-"] or (len(args.genes) == 1 and args.genes[0] == "-"):
        genes = [ln.strip() for ln in sys.stdin if ln.strip()]
    else:
        genes = [g.strip().upper() for g in args.genes if g.strip()]

    if not genes:
        print("No genes provided", file=sys.stderr)
        return 1

    for g in genes:
        print(f"Fetching {g}...", flush=True)
        d = get_or_fetch_gene_diseases(
            g,
            args.db,
            api_key,
            model=args.model,
            allow_gemini=True,
        )
        if d:
            print(f"  OK: {d[0].get('name')} ({d[0].get('inheritance')})")
        else:
            print(f"  FAILED: {g}", file=sys.stderr)

    print(f"Done. Database: {args.db}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
