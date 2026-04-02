"""
Configurable WES / exome interpretation panels (order-time selection).

- **Bundled** catalog: ``data/wes_panels.json`` (or ``WES_PANELS_JSON``).
- **Custom** (portal-built): ``data/wes_panels_custom.json`` (or ``WES_PANELS_CUSTOM_JSON``); merged; custom wins on id clash.

Panels should define ``interpretation_genes`` (HGNC symbols) for report scope; the carrier plugin
narrows ``result.json`` to that set **after** annotation / ACMG (``panel_filter_after_analysis``).
Optional ``disease_bed`` / ``backbone_bed`` are for capture regions and ``disease_bed_gene`` tagging;
coordinate pre-filter from disease BED is **off** by default (``restrict_to_disease_bed``). Order params may add
``interpretation_genes_extra`` for reflex testing without re-running the pipeline.
"""

from __future__ import annotations

import gzip
import json
import logging
import os
import re
from typing import Any, Dict, List, Optional, Set, Tuple

from ..config import settings

logger = logging.getLogger(__name__)

_ID_RE = re.compile(r"^[a-z][a-z0-9_-]{1,62}$")


def _project_root_dir() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _default_catalog_path() -> str:
    env = (settings.wes_panels_json or "").strip()
    if env:
        return env if os.path.isabs(env) else os.path.normpath(os.path.join(_project_root_dir(), env))
    return os.path.join(_project_root_dir(), "data", "wes_panels.json")


def _should_use_data_wes_panels_fallback(abs_default: str) -> bool:
    """
    Docker often mounts ``/app`` or ``/app/data`` read-only; default repo paths then fail with
    Errno 30. When ``/data`` is mounted rw (typical compose), use that instead.
    """
    if not abs_default.startswith("/app/"):
        return False
    parent = os.path.dirname(abs_default)
    try:
        if os.path.isdir(parent) and os.access(parent, os.W_OK):
            return False
    except OSError:
        pass
    try:
        return bool(os.path.isdir("/data") and os.access("/data", os.W_OK))
    except OSError:
        return False


def _custom_catalog_path() -> str:
    env = (settings.wes_panels_custom_json or "").strip()
    if env:
        return env if os.path.isabs(env) else os.path.normpath(os.path.join(_project_root_dir(), env))
    default = os.path.join(_project_root_dir(), "data", "wes_panels_custom.json")
    if _should_use_data_wes_panels_fallback(os.path.normpath(default)):
        return "/data/wes_panels/wes_panels_custom.json"
    return default


def _generated_dir() -> str:
    env = (settings.wes_panels_generated_dir or "").strip()
    if env:
        return env if os.path.isabs(env) else os.path.normpath(os.path.join(_project_root_dir(), env))
    default = os.path.join(_project_root_dir(), "data", "bed", "wes_panels_generated")
    if _should_use_data_wes_panels_fallback(os.path.normpath(default)):
        return "/data/wes_panels/generated"
    return default


def _resolve_repo_path(p: str) -> str:
    p = (p or "").strip()
    if not p:
        return ""
    if os.path.isabs(p):
        return os.path.normpath(p)
    return os.path.normpath(os.path.join(_project_root_dir(), p))


def _load_panels_array(path: str) -> List[Dict[str, Any]]:
    if not os.path.isfile(path):
        return []
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        raw = data.get("panels") if isinstance(data, dict) else None
        if not isinstance(raw, list):
            return []
        # Coerce id to str: JSON may have numeric id; (n or "").strip() would raise AttributeError.
        return [p for p in raw if isinstance(p, dict) and str(p.get("id") or "").strip()]
    except Exception as e:
        logger.error("[wes_panels] Cannot read %s: %s", path, e)
        return []


def _bundled_ids() -> Set[str]:
    return {str(p.get("id", "")).strip() for p in _load_panels_array(_default_catalog_path()) if p.get("id")}


def load_wes_panels_raw() -> Dict[str, Any]:
    """Bundled file only (for admin / version)."""
    path = _default_catalog_path()
    if not os.path.isfile(path):
        logger.warning("[wes_panels] Catalog not found: %s", path)
        return {"version": 1, "panels": []}
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        logger.error("[wes_panels] Cannot read %s: %s", path, e)
        return {"version": 1, "panels": []}


def list_panels() -> List[Dict[str, Any]]:
    """Merged panels: bundled first, then custom (custom overrides same id)."""
    by_id: Dict[str, Dict[str, Any]] = {}
    for p in _load_panels_array(_default_catalog_path()):
        pid = str(p.get("id", "")).strip()
        if not pid:
            continue
        entry = dict(p)
        entry["_source"] = "bundled"
        by_id[pid] = entry
    for p in _load_panels_array(_custom_catalog_path()):
        pid = str(p.get("id", "")).strip()
        if not pid:
            continue
        entry = dict(p)
        entry["_source"] = "custom"
        by_id[pid] = entry
    return list(by_id.values())


def get_panel_by_id(panel_id: Optional[str]) -> Optional[Dict[str, Any]]:
    pid = (panel_id or "").strip()
    if not pid:
        return None
    for p in list_panels():
        if str(p.get("id", "")).strip() == pid:
            return p
    return None


_GENE_SYMBOL_LIKE = re.compile(r"^[A-Z][A-Z0-9-]{0,62}$")
# Typical HGNC symbols are short; BEDs may label syndromes / regions in column 4 (e.g. WOLF-HIRSCHHORN).
_MAX_HGNC_LIKE_LEN = 15


def _hgnc_like_gene_symbol(name: str) -> bool:
    """
    True if the BED 4th-column name looks like an HGNC gene symbol (not cytoband like 1p36 or 2p25.3).
    Excludes long syndrome-style labels and multi-hyphen region names.
    """
    s = (name or "").strip().upper()
    if not s:
        return False
    if "." in s:
        return False
    if len(s) > _MAX_HGNC_LIKE_LEN:
        return False
    if s.count("-") > 1:
        return False
    return bool(_GENE_SYMBOL_LIKE.match(s))


def load_gene_symbols_from_text_file(abs_path: str) -> List[str]:
    """One gene per line; # comments; comma/semicolon separated lines allowed."""
    if not abs_path or not os.path.isfile(abs_path):
        return []
    out: List[str] = []
    try:
        with open(abs_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = (line or "").split("#", 1)[0].strip()
                if not line:
                    continue
                for tok in re.split(r"[\s,;]+", line):
                    t = _clean_gene_token(tok)
                    if t:
                        out.append(t.upper())
    except OSError as e:
        logger.warning("[wes_panels] Cannot read interpretation_genes_file %s: %s", abs_path, e)
        return []
    return list(dict.fromkeys(out))


def load_gene_symbols_from_bed_column4(abs_path: str) -> List[str]:
    """
    Unique HGNC-like symbols from BED column 4 (excludes cytobands and region labels).
    Use when ``interpretation_genes_from_disease_bed`` is true: one technical BED, report scope = genes.
    """
    if not abs_path or not os.path.isfile(abs_path):
        return []
    seen: set = set()
    out: List[str] = []
    try:
        opener = gzip.open if abs_path.endswith(".gz") else open
        with opener(abs_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                raw = (parts[3] or "").strip()
                if not _hgnc_like_gene_symbol(raw):
                    continue
                u = raw.upper()
                if u not in seen:
                    seen.add(u)
                    out.append(u)
    except OSError as e:
        logger.warning("[wes_panels] Cannot read disease BED for gene list %s: %s", abs_path, e)
        return []
    return out


def resolve_panel_interpretation_genes(panel: Dict[str, Any]) -> List[str]:
    """
    Merged, deduped HGNC symbols for the panel (for post-annotation filtering).

    Sources (in order):
    - ``interpretation_genes`` array in JSON
    - ``interpretation_genes_file`` (repo-relative path, one gene per line)
    - ``interpretation_genes_from_disease_bed``: extract HGNC-like symbols from ``disease_bed`` column 4
    """
    seen: set = set()
    ordered: List[str] = []

    def add_many(xs: List[str]) -> None:
        for x in xs:
            u = str(x).strip().upper()
            if not u:
                continue
            if u not in seen:
                seen.add(u)
                ordered.append(u)

    raw_list = panel.get("interpretation_genes")
    if isinstance(raw_list, list):
        add_many([str(x).strip() for x in raw_list if str(x).strip()])

    gf = str(panel.get("interpretation_genes_file") or "").strip()
    if gf:
        abs_gf = _resolve_repo_path(gf)
        add_many(load_gene_symbols_from_text_file(abs_gf))

    if bool(panel.get("interpretation_genes_from_disease_bed")):
        rel = str(panel.get("disease_bed") or "").strip()
        if rel:
            abs_bed = _resolve_repo_path(rel)
            add_many(load_gene_symbols_from_bed_column4(abs_bed))
        else:
            logger.warning(
                "[wes_panels] Panel %s has interpretation_genes_from_disease_bed but no disease_bed",
                panel.get("id"),
            )

    return ordered


def approximate_gene_count_from_bed(abs_path: str) -> Optional[int]:
    """Rough count: unique non-empty 4th column tokens, else non-comment line count."""
    if not abs_path or not os.path.isfile(abs_path):
        return None
    genes: set = set()
    lines = 0
    try:
        with open(abs_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                lines += 1
                parts = line.split("\t")
                if len(parts) >= 4:
                    g = (parts[3] or "").strip()
                    if g:
                        genes.add(g)
        if genes:
            return len(genes)
        return lines if lines else None
    except OSError:
        return None


def interpretation_gene_set_for_job(job: Any) -> Set[str]:
    """
    Genes to retain in the interpretation report when post-filtering is active.
    Merges panel list (``panel_interpretation_genes``) with order extras
    (``interpretation_genes_extra`` and optional ``carrier.interpretation_genes_extra``).
    """
    out: Set[str] = set()
    params = getattr(job, "params", None) or {}
    raw = params.get("panel_interpretation_genes")
    if isinstance(raw, list):
        out.update(str(x).strip().upper() for x in raw if str(x).strip())
    elif isinstance(raw, str) and raw.strip():
        out.update(g.strip().upper() for g in raw.split(",") if g.strip())

    def _merge_extra(txt: Optional[str]) -> None:
        if not txt or not str(txt).strip():
            return
        out.update(g.strip().upper() for g in str(txt).split(",") if g.strip())

    _merge_extra(params.get("interpretation_genes_extra"))
    carrier = params.get("carrier")
    if isinstance(carrier, dict):
        _merge_extra(carrier.get("interpretation_genes_extra"))
    return out


def should_apply_interpretation_post_filter(job: Any) -> bool:
    """True when result.json should list only variants whose gene is in the interpretation set."""
    params = getattr(job, "params", None) or {}
    if params.get("panel_filter_after_analysis") is False:
        return False
    genes = interpretation_gene_set_for_job(job)
    if not genes:
        return False
    if params.get("panel_filter_after_analysis") is True:
        return True
    ex = (params.get("interpretation_genes_extra") or "").strip()
    if ex:
        return True
    carrier = params.get("carrier")
    if isinstance(carrier, dict) and (carrier.get("interpretation_genes_extra") or "").strip():
        return True
    return False


def panels_for_api_response() -> List[Dict[str, Any]]:
    """Sanitized list for the portal (paths resolved to absolute for display)."""
    out: List[Dict[str, Any]] = []
    for p in list_panels():
        pid = str(p.get("id", "")).strip()
        disease = _resolve_repo_path(str(p.get("disease_bed") or ""))
        backbone = _resolve_repo_path(str(p.get("backbone_bed") or ""))
        gc = p.get("gene_count")
        ig = p.get("interpretation_genes")
        resolved = resolve_panel_interpretation_genes(p)
        if gc is None and isinstance(ig, list) and ig:
            gc = len([x for x in ig if str(x).strip()])
        if gc is None and resolved:
            gc = len(resolved)
        if gc is None and disease and os.path.isfile(disease):
            gc = approximate_gene_count_from_bed(disease)
        src = p.get("_source") or "bundled"
        item: Dict[str, Any] = {
            "id": pid,
            "label": p.get("label") or pid,
            "category": p.get("category") or "other",
            "description": p.get("description"),
            "gene_count": gc,
            "disease_bed_resolved": disease if os.path.isfile(disease) else None,
            "backbone_bed_resolved": backbone if backbone and os.path.isfile(backbone) else None,
            "source": src if src in ("bundled", "custom") else "bundled",
        }
        out.append(item)
    out.sort(key=lambda x: (x.get("category") or "", x.get("label") or ""))
    return out


def apply_wes_panel_to_job_params(job: Any) -> None:
    """
    Merge catalog BED paths into ``job.params`` when ``wes_panel_id`` is set.
    Explicit ``disease_bed`` / ``backbone_bed`` in params always win.

    If the panel defines ``interpretation_genes``, those symbols are stored in
    ``panel_interpretation_genes`` and ``panel_filter_after_analysis`` is set so the
    carrier plugin can narrow ``result.json`` *after* parse/ACMG. The panel's
    ``disease_bed`` is applied only when ``narrow_with_panel_bed`` is true (legacy
    panels without ``interpretation_genes`` behave as before: apply disease BED when set).
    """
    if not getattr(job, "params", None):
        job.params = {}
    assert job.params is not None

    pid = job.params.get("wes_panel_id")
    if not pid and isinstance(job.params.get("carrier"), dict):
        pid = job.params["carrier"].get("wes_panel_id")
    pid = (pid or "").strip()
    if not pid:
        return

    panel = get_panel_by_id(pid)
    if not panel:
        logger.warning("[wes_panels] Unknown wes_panel_id=%r — ignoring", pid)
        return

    interp_list = resolve_panel_interpretation_genes(panel)
    narrow_from_panel_bed = True
    if interp_list:
        job.params["panel_interpretation_genes"] = interp_list
        job.params["panel_filter_after_analysis"] = True
        narrow_from_panel_bed = bool(panel.get("narrow_with_panel_bed", False))
        logger.info(
            "[wes_panels] Panel %s: interpretation_genes=%d, narrow_with_panel_bed=%s",
            pid,
            len(interp_list),
            narrow_from_panel_bed,
        )
    else:
        job.params.pop("panel_interpretation_genes", None)
        # Leave panel_filter_after_analysis as-is if set only by order (extras-only path).

    if not (job.params.get("disease_bed") or "").strip():
        rel = str(panel.get("disease_bed") or "").strip()
        if rel and narrow_from_panel_bed:
            abs_p = _resolve_repo_path(rel)
            if os.path.isfile(abs_p):
                job.params["disease_bed"] = abs_p
                logger.info("[wes_panels] Applied disease_bed from panel %s → %s", pid, abs_p)
            else:
                logger.warning("[wes_panels] Panel %s disease_bed not found: %s", pid, abs_p)
        elif rel and not narrow_from_panel_bed:
            logger.info(
                "[wes_panels] Skipping panel %s disease_bed (interpretation_genes + narrow_with_panel_bed=false)",
                pid,
            )

    if not (job.params.get("backbone_bed") or "").strip():
        rel = str(panel.get("backbone_bed") or "").strip()
        if rel:
            abs_p = _resolve_repo_path(rel)
            if os.path.isfile(abs_p):
                job.params["backbone_bed"] = abs_p
                logger.info("[wes_panels] Applied backbone_bed from panel %s → %s", pid, abs_p)
            else:
                logger.warning("[wes_panels] Panel %s backbone_bed not found: %s", pid, abs_p)


def normalize_panel_id(raw: str) -> str:
    s = (raw or "").strip().lower()
    s = re.sub(r"[^a-z0-9_-]+", "_", s)
    s = s.strip("_")
    if not _ID_RE.match(s):
        raise ValueError(
            "Panel id must be 2–64 chars: start with a letter, then letters, digits, _ or - "
            f"(got {raw!r} → {s!r})"
        )
    return s


def _clean_gene_token(raw: str) -> str:
    """
    Normalize a pasted gene symbol: strip BOM/whitespace and trailing list punctuation
    (e.g. ``ZNF469.`` at end of a sentence).
    """
    s = (raw or "").strip().replace("\ufeff", "").replace("\u00a0", " ")
    s = s.strip(".,;:\t\r\n")
    s = s.rstrip(".")
    return s.strip()


def split_genes_from_text(text: Optional[str]) -> List[str]:
    if not text or not str(text).strip():
        return []
    raw = str(text).replace("\ufeff", "").replace("\u00a0", " ")
    # Fullwidth comma/semicolon (common in pasted Office/docs text)
    raw = raw.replace("\uff0c", ",").replace("\uff1b", ";")
    parts = re.split(r"[\s,;]+", raw)
    out: List[str] = []
    for p in parts:
        t = _clean_gene_token(p)
        if not t:
            continue
        if len(t) > 64:
            raise ValueError(
                "Paste uses commas or line breaks between gene symbols. "
                f"One token is too long to be a symbol (starts with {t[:40]!r}…)."
            )
        out.append(t)
    return out


def _coerce_gene_list_from_ambiguous_disease_field(raw: Optional[str]) -> Optional[List[str]]:
    """
    If the user pasted a gene list into the "disease BED path" field (or Docker prepended /app/…),
    detect that and return parsed symbols instead of treating the blob as a path.
    """
    s = (raw or "").strip()
    if not s or len(s) < 40:
        return None
    abs_try = os.path.normpath(s) if os.path.isabs(s) else _resolve_repo_path(s)
    if os.path.isfile(abs_try) or os.path.isfile(s):
        return None
    if s.lower().endswith(".bed") and "\n" not in s and s.count(",") <= 2 and len(s) < 800:
        return None
    if "," not in s and "\n" not in s and "\r" not in s:
        return None

    _junk = frozenset(
        {"app", "data", "bed", "src", "usr", "var", "home", "dev", "twist", "opt", "tmp", "root"}
    )
    pat = re.compile(r"^[A-Za-z][A-Za-z0-9\-]*$")
    chunks = re.split(r"[\s,;]+", s)
    seen: Set[str] = set()
    cleaned: List[str] = []
    for chunk in chunks:
        if not chunk:
            continue
        t = re.sub(r"^.*[/\\\\]", "", chunk.strip())
        t = _clean_gene_token(t)
        if not pat.match(t) or len(t) > 40:
            continue
        if t.lower() in _junk:
            continue
        u = t.upper()
        if u in seen:
            continue
        seen.add(u)
        cleaned.append(t)
    if len(cleaned) >= 5:
        return cleaned
    return None


def resolve_gene_source_bed(explicit: Optional[str]) -> str:
    """First usable path: explicit (portal), then WES_PANEL_GENE_SOURCE_BED, then GENE_BED."""
    candidates = [
        (explicit or "").strip(),
        (settings.wes_panel_gene_source_bed or "").strip(),
        (settings.gene_bed or "").strip(),
    ]
    for c in candidates:
        if not c:
            continue
        abs_p = os.path.normpath(c) if os.path.isabs(c) else _resolve_repo_path(c)
        if os.path.isfile(abs_p):
            return abs_p
    return ""


def build_disease_bed_from_genes(
    genes: List[str],
    gene_bed_path: str,
    out_abs_path: str,
) -> Tuple[int, List[str], List[str]]:
    """
    Subset ``gene_bed_path`` (chrom, start, end, name) to rows whose 4th column matches
    requested gene symbols (case-insensitive).

    Returns (lines_written, matched_genes_sorted, missing_genes_sorted).
    """
    wanted: Set[str] = {g.strip().upper() for g in genes if g and g.strip()}
    if not wanted:
        raise ValueError("No genes provided")
    if not gene_bed_path or not os.path.isfile(gene_bed_path):
        raise ValueError(
            f"Source BED for gene intervals not found or not configured: {gene_bed_path!r}. "
            "Set GENE_BED or WES_PANEL_GENE_SOURCE_BED (e.g. Twist capture BED with gene symbol in column 4)."
        )

    matched_lines: List[str] = []
    matched_names: Set[str] = set()

    with open(gene_bed_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            raw = line.rstrip("\n\r")
            if not raw.strip() or raw.startswith("#") or raw.startswith("track"):
                continue
            parts = raw.split("\t")
            if len(parts) < 4:
                continue
            name = (parts[3] or "").strip()
            if name.upper() in wanted:
                matched_lines.append(raw)
                matched_names.add(name.upper())

    missing = sorted(wanted - matched_names)
    matched_sorted = sorted(matched_names)

    os.makedirs(os.path.dirname(out_abs_path) or ".", exist_ok=True)
    header = (
        f"# wes_panels generated — {len(matched_lines)} intervals, {len(matched_sorted)} genes\n"
        f"# source_gene_bed={gene_bed_path}\n"
    )
    with open(out_abs_path, "w", encoding="utf-8") as out:
        out.write(header)
        for ln in matched_lines:
            out.write(ln + "\n")

    return len(matched_lines), [m for m in matched_sorted], missing


def _rel_path_for_json(abs_path: str) -> str:
    """Store repo-relative path when under project root for portability."""
    root = _project_root_dir()
    try:
        ap = os.path.normpath(os.path.abspath(abs_path))
        rp = os.path.normpath(root)
        if ap.startswith(rp + os.sep):
            return os.path.relpath(ap, root).replace("\\", "/")
    except OSError:
        pass
    return abs_path


def save_custom_panel(
    panel_id: str,
    label: str,
    category: str,
    description: Optional[str],
    backbone_bed: Optional[str],
    disease_bed: Optional[str],
    genes: Optional[List[str]],
    genes_text: Optional[str],
    gene_source_bed: Optional[str] = None,
    skip_generated_bed: bool = False,
) -> Dict[str, Any]:
    """
    Write or replace one panel in the custom JSON file.

    Gene list mode: either build a subset BED from a source interval BED (when
    ``skip_generated_bed`` is false and a source BED resolves), or store
    ``interpretation_genes`` only for post-analysis filtering (no generated file).
    Alternatively provide an existing ``disease_bed`` path (BED-file mode).
    """
    pid = normalize_panel_id(panel_id)
    if pid in _bundled_ids():
        raise ValueError(f"Id {pid!r} is reserved for a built-in panel — choose another id")

    seen: Set[str] = set()
    cleaned: List[str] = []
    for g in (genes or []) + split_genes_from_text(genes_text):
        t = _clean_gene_token(str(g))
        if not t:
            continue
        u = t.upper()
        if u in seen:
            continue
        seen.add(u)
        cleaned.append(u)

    db_path_in: Optional[str] = (disease_bed or "").strip() or None

    # User may have pasted the gene list into the BED path field — treat as genes, not a path.
    if not cleaned and db_path_in:
        coerced = _coerce_gene_list_from_ambiguous_disease_field(db_path_in)
        if coerced:
            cleaned = coerced
            db_path_in = None

    rel_disease: Optional[str] = None
    matched: Optional[List[str]] = None
    warn_missing: Optional[List[str]] = None
    source_rel: Optional[str] = None

    if cleaned:
        gb = "" if skip_generated_bed else resolve_gene_source_bed(gene_source_bed)
        if not skip_generated_bed and gb:
            out_abs = os.path.join(_generated_dir(), f"{pid}.bed")
            n_lines, matched, missing = build_disease_bed_from_genes(cleaned, gb, out_abs)
            if n_lines == 0:
                raise ValueError(
                    "No intervals matched your gene list in the source BED. "
                    "Check gene symbols match column 4 of that BED, or use a BED that lists HGNC symbols."
                )
            rel_disease = _rel_path_for_json(out_abs)
            warn_missing = missing
            source_rel = _rel_path_for_json(gb)
        else:
            matched = []
            warn_missing = []
            if not skip_generated_bed and not gb:
                raise ValueError(
                    "Building a panel BED from a gene list requires a source interval BED "
                    "(gene symbol in column 4). Set GENE_BED or WES_PANEL_GENE_SOURCE_BED, "
                    "enter “Source BED” in the portal, or enable “interpretation gene list only” "
                    "to skip generating a file."
                )
    elif db_path_in:
        abs_db = _resolve_repo_path(db_path_in) if not os.path.isabs(db_path_in) else os.path.normpath(db_path_in)
        if not os.path.isfile(abs_db):
            raise ValueError(f"disease_bed file not found: {abs_db}")
        rel_disease = _rel_path_for_json(abs_db)
        matched = []
        warn_missing = []
        source_rel = None
    else:
        raise ValueError(
            "No gene symbols were parsed. Use commas or line breaks between symbols "
            "(check for a stray period after the last gene, or paste from plain text)."
        )

    bb_rel: Optional[str] = None
    if (backbone_bed or "").strip():
        abs_bb = _resolve_repo_path(backbone_bed) if not os.path.isabs(backbone_bed) else os.path.normpath(backbone_bed)
        if not os.path.isfile(abs_bb):
            raise ValueError(f"backbone_bed file not found: {abs_bb}")
        bb_rel = _rel_path_for_json(abs_bb)

    panel_obj: Dict[str, Any] = {
        "id": pid,
        "label": label.strip(),
        "category": (category or "other").strip() or "other",
        "description": (description or "").strip() or None,
    }
    if rel_disease:
        panel_obj["disease_bed"] = rel_disease
    if cleaned:
        panel_obj["interpretation_genes"] = [str(g).strip() for g in cleaned]
        panel_obj["narrow_with_panel_bed"] = False
        panel_obj["gene_count"] = len(cleaned)
    if bb_rel:
        panel_obj["backbone_bed"] = bb_rel
    if panel_obj.get("description") is None:
        panel_obj.pop("description", None)
    if cleaned and source_rel:
        panel_obj["gene_interval_source_bed"] = source_rel

    custom_path = _custom_catalog_path()
    os.makedirs(os.path.dirname(custom_path) or ".", exist_ok=True)

    existing = _load_panels_array(custom_path)
    others = [p for p in existing if str(p.get("id", "")).strip() != pid]
    others.append(panel_obj)
    others.sort(key=lambda x: str(x.get("id", "")))

    payload = {"version": 1, "panels": others}
    with open(custom_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
        f.write("\n")

    logger.info("[wes_panels] Saved custom panel %s → %s", pid, custom_path)

    out: Dict[str, Any] = {
        "id": pid,
        "disease_bed": rel_disease,
        "genes_matched": matched if cleaned else None,
        "genes_missing_from_bed": warn_missing if cleaned else None,
        "interpretation_genes_only": bool(cleaned) and rel_disease is None,
    }
    if cleaned and source_rel:
        out["gene_interval_source_bed"] = source_rel
    return out


def delete_custom_panel(panel_id: str) -> bool:
    """Remove from custom JSON only. Returns False if not found or id is bundled-only."""
    pid = normalize_panel_id(panel_id)
    if pid in _bundled_ids():
        raise ValueError("Cannot delete a built-in panel id")

    custom_path = _custom_catalog_path()
    existing = _load_panels_array(custom_path)
    found = None
    rest: List[Dict[str, Any]] = []
    for p in existing:
        if str(p.get("id", "")).strip() == pid:
            found = p
        else:
            rest.append(p)
    if not found:
        return False

    gen = os.path.join(_generated_dir(), f"{pid}.bed")
    if os.path.isfile(gen):
        try:
            os.remove(gen)
        except OSError as e:
            logger.warning("[wes_panels] Could not remove generated BED %s: %s", gen, e)

    payload = {"version": 1, "panels": rest}
    with open(custom_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
        f.write("\n")
    logger.info("[wes_panels] Deleted custom panel %s", pid)
    return True
