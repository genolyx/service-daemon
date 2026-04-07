#!/usr/bin/env python3
"""
Compare two carrier review ``result.json`` files (``variants[]`` from portal pipeline).

Usage:
  python scripts/compare_carrier_result_json.py \\
    /path/to/acmg_result.json /path/to/invitae_result.json \\
    --genes-acmg data/gene_lists/acmg_genes.txt \\
    --genes-invitae data/gene_lists/invitae_genes.txt

Optional gene lists:
  One HGNC symbol per line (``#`` comments ok). Used to label variants that appear
  in only one run: **only on A panel**, **only on B panel**, **both panels**, **neither**.

If you omit ``--genes-*``, the script still prints variant/gene differences; you can
classify genes manually or load lists from ``resolve_panel_interpretation_genes`` in
the daemon environment.

Example (panel lists from the running app):
  cd /path/to/service-daemon && PYTHONPATH=. python -c "
  from app.services.wes_panels import get_panel_by_id, resolve_panel_interpretation_genes
  for pid in ('invitae_302',):
      g = resolve_panel_interpretation_genes(get_panel_by_id(pid) or {})
      open(f'/tmp/{pid}_genes.txt','w').write('\\n'.join(g))
  "
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple


def _norm_chrom(c: Any) -> str:
    s = str(c or "").strip()
    if not s:
        return ""
    if s.startswith("chr"):
        return s
    if s in ("X", "Y", "M", "MT"):
        return "chr" + s
    return "chr" + s if s.isdigit() or (len(s) <= 2 and s.startswith(("X", "Y"))) else s


def variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    return "|".join(
        [
            _norm_chrom(chrom),
            str(int(pos)) if str(pos).isdigit() else str(pos),
            str(ref or "").upper(),
            str(alt or "").upper(),
        ]
    )


def load_genes(path: Optional[str]) -> Set[str]:
    if not path or not os.path.isfile(path):
        return set()
    out: Set[str] = set()
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            for tok in re.split(r"[\s,;]+", line):
                t = tok.strip().upper()
                if t:
                    out.add(t)
    return out


def load_variants(path: str) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    with open(path, encoding="utf-8", errors="replace") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise SystemExit(f"Expected JSON object in {path}")
    raw = data.get("variants")
    if not isinstance(raw, list):
        raw = []
    return raw, data


def index_by_key(variants: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for v in variants:
        if not isinstance(v, dict):
            continue
        k = variant_key(v.get("chrom"), v.get("pos"), v.get("ref"), v.get("alt"))
        out[k] = v
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare two carrier result.json variants[]")
    ap.add_argument("a", help="result.json path (e.g. ACMG carrier run)")
    ap.add_argument("b", help="result.json path (e.g. Invitae_302 run)")
    ap.add_argument("--label-a", default="A", help="label for first file")
    ap.add_argument("--label-b", default="B", help="label for second file")
    ap.add_argument(
        "--genes-a",
        metavar="FILE",
        help="HGNC gene list for panel A (e.g. ACMG interpretation genes)",
    )
    ap.add_argument(
        "--genes-b",
        metavar="FILE",
        help="HGNC gene list for panel B (e.g. Invitae_302 interpretation genes)",
    )
    args = ap.parse_args()

    va, ja = load_variants(args.a)
    vb, jb = load_variants(args.b)

    ga = load_genes(args.genes_a)
    gb = load_genes(args.genes_b)

    ia = index_by_key(va)
    ib = index_by_key(vb)
    ka = set(ia.keys())
    kb = set(ib.keys())

    def _gene(v: Dict[str, Any]) -> str:
        return str(v.get("gene") or "").strip().upper()

    print("=== filter_summary (if present) ===")
    for label, j in ((args.label_a, ja), (args.label_b, jb)):
        fs = j.get("filter_summary") if isinstance(j.get("filter_summary"), dict) else {}
        print(f"\n{label}:")
        for k in (
            "wes_panel_id",
            "wes_panel_label",
            "interpretation_post_filter_genes",
            "interpretation_post_filter_applied",
            "backbone_bed",
            "disease_bed",
            "restrict_to_disease_bed",
            "require_protein_altering",
        ):
            if k in fs:
                print(f"  {k}: {fs.get(k)}")

    print("\n=== counts ===")
    print(f"{args.label_a} variants: {len(va)}")
    print(f"{args.label_b} variants: {len(vb)}")
    print(f"Unique keys {args.label_a}: {len(ka)}")
    print(f"Unique keys {args.label_b}: {len(kb)}")
    print(f"Intersection (same chrom/pos/ref/alt): {len(ka & kb)}")
    only_a = ka - kb
    only_b = kb - ka
    print(f"Only in {args.label_a}: {len(only_a)}")
    print(f"Only in {args.label_b}: {len(only_b)}")

    def classify_gene(gene: str) -> str:
        g = gene.strip().upper()
        if not g:
            return "empty_gene"
        in_a = g in ga if ga else None
        in_b = g in gb if gb else None
        if in_a is None and in_b is None:
            return "no_gene_lists"
        if in_a and in_b:
            return "both_panels"
        if in_a:
            return f"only_{args.label_a}_panel"
        if in_b:
            return f"only_{args.label_b}_panel"
        return "neither_panel"

    def report_side(prefix: str, keys: Set[str], src: Dict[str, Dict[str, Any]]) -> None:
        by_gene: Dict[str, int] = {}
        by_class: Dict[str, int] = {}
        for k in sorted(keys):
            v = src.get(k) or {}
            g = _gene(v)
            by_gene[g] = by_gene.get(g, 0) + 1
            c = classify_gene(g)
            by_class[c] = by_class.get(c, 0) + 1
        print(f"\n=== {prefix} ({len(keys)} variants) ===")
        if by_class:
            print("By panel membership (gene symbol):")
            for c in sorted(by_class.keys(), key=lambda x: (-by_class[x], x)):
                print(f"  {c}: {by_class[c]}")
        print("By gene (count):")
        for g in sorted(by_gene.keys(), key=lambda x: (-by_gene[x], x)):
            print(f"  {g or '(no gene)'}: {by_gene[g]}")

    report_side(f"Only in {args.label_a}", only_a, ia)
    report_side(f"Only in {args.label_b}", only_b, ib)

    # Intersection: same variant, but annotation might differ — optional diff
    both = ka & kb
    if both:
        diff_fields = []
        for k in sorted(both)[:5000]:
            a, b = ia[k], ib[k]
            for fld in [
                "gene",
                "hgvsc",
                "hgvsp",
                "effect",
                "clinvar_sig_primary",
                "gnomad_af",
            ]:
                if (a.get(fld) or "") != (b.get(fld) or ""):
                    diff_fields.append((k, fld, a.get(fld), b.get(fld)))
        if diff_fields:
            print(f"\n=== same variant key, different field (first 20) ===")
            for row in diff_fields[:20]:
                print(f"  {row[0]} | {row[1]}: {args.label_a}={row[2]!r} {args.label_b}={row[3]!r}")
        else:
            print(f"\n=== same variant key: no gene/hgv/clinvar/gnomad diffs on {len(both)} pairs ===")

    if not ga and not gb:
        print(
            "\nTip: pass --genes-a and --genes-b to label only-in-A / only-in-B variants "
            "by panel membership.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
