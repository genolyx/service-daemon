"""
Carrier Screening Service Plugin

ServicePlugin 인터페이스 구현체.
전체 워크플로우를 조율합니다:

    1. prepare_inputs     : FASTQ 파일 확인/다운로드, 디렉토리 구조 생성
    2. get_pipeline_command: Nextflow 파이프라인 실행 명령 생성
    3. check_completion   : 파이프라인 완료 확인 (VCF 존재 여부)
    4. process_results    : VCF → BED 필터 → Annotation → ACMG → result.json 생성
    5. get_output_files   : Portal에 업로드할 파일 목록 반환
    6. generate_report    : 리뷰어 확정 후 report.json + 다국어 PDF 생성
"""

import os
import re
import glob
import asyncio
import json
import logging
import shlex
import shutil
from typing import Any, Dict, List, Optional, Set, Tuple

from ..base import ServicePlugin
from ..wes_panels import (
    interpretation_gene_set_for_job,
    pgx_portal_gene_allowlist_for_job,
    should_apply_interpretation_post_filter,
)
from ...config import normalize_legacy_carrier_container_path, settings
from ...models import Job, OutputFile
from .layout_norm import (
    carrier_run_analysis_work_arg,
    carrier_sequencing_folder,
)

logger = logging.getLogger(__name__)

_GX_EXOME_LAYOUT_PREFIX = "/data/gx-exome"


def _resolve_nextflow_executable() -> str:
    """
    Avoid exit 127 when NEXTFLOW_EXECUTABLE=nextflow but PATH is minimal in a subprocess.
    Prefer a real file path (image default: /usr/local/bin/nextflow).
    """
    exe = (settings.nextflow_executable or "").strip() or "nextflow"
    if os.path.isabs(exe) and os.path.isfile(exe):
        return exe
    if not os.path.isabs(exe) and os.path.isfile(exe):
        return os.path.abspath(exe)
    found = shutil.which(exe)
    if found:
        return found
    for cand in ("/usr/local/bin/nextflow", "/usr/bin/nextflow"):
        if os.path.isfile(cand):
            return cand
    raise RuntimeError(
        "Nextflow executable not found. Set NEXTFLOW_EXECUTABLE to the full path "
        "(e.g. NEXTFLOW_EXECUTABLE=/usr/local/bin/nextflow in .env.docker) or install Nextflow on PATH."
    )


def _bash_executable() -> str:
    """bash for run_analysis.sh; fail fast with a clear error instead of exit 127."""
    w = shutil.which("bash")
    if w and os.path.isfile(w):
        return w
    for cand in ("/bin/bash", "/usr/bin/bash"):
        if os.path.isfile(cand):
            return cand
    raise RuntimeError(
        "bash not found (expected /bin/bash or /usr/bin/bash). "
        "The service-daemon image should include bash; check PATH or base image."
    )


def carrier_report_generation_paths(job: Job) -> Tuple[str, str]:
    """
    (output_dir, analysis_dir) for Generate Report: report.json/PDF and QC extraction.

    When settings.carrier_screening_report_output_root is set, uses
    <root>/output|analysis/<work_dir>/<sample_name> under that tree (e.g. /home/sam/Carrier_result).
    Otherwise uses this order's work_dir + sample_name under carrier_screening_work_root.

    Do **not** use job.output_dir here: ``align_carrier_job_dirs_from_main_vcf`` may point
    ``job.output_dir`` at a *prior* pipeline tree when reusing the same VCF; portal
    ``result.json`` must still live under the **current** order's folder.
    """
    root = (settings.carrier_screening_report_output_root or "").strip()
    wk = str(job.work_dir).strip() or "00"
    smp = str(job.sample_name).strip()
    if root:
        return (
            os.path.join(root, "output", wk, smp),
            os.path.join(root, "analysis", wk, smp),
        )
    work_root = settings.carrier_screening_work_root
    canonical_out = os.path.join(work_root, "output", wk, smp)
    canonical_an = os.path.join(work_root, "analysis", wk, smp)
    return canonical_out, canonical_an


def carrier_report_output_dir(job: Job) -> str:
    """Directory where report.json and Report_*.pdf are written for this order."""
    return carrier_report_generation_paths(job)[0]


def _extra_result_json_paths_for_carrier_report(job: Job) -> List[str]:
    """
    When report output root differs from pipeline ``job.output_dir``, ``result.json`` may only
    exist under the work tree. Try that path after the report-directory result.json.
    """
    out: List[str] = []
    primary = os.path.normpath(
        os.path.abspath(os.path.join(carrier_report_output_dir(job), "result.json"))
    )
    if job.output_dir:
        alt = os.path.join(job.output_dir, "result.json")
        if os.path.normpath(os.path.abspath(alt)) != primary:
            out.append(alt)
    return out


def _filter_variants_by_interpretation_genes(
    annotated_variants: List[Dict[str, Any]],
    acmg_results: List[Dict[str, Any]],
    gene_set: Set[str],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Keep variants whose annotated gene symbol is in ``gene_set`` (HGNC uppercase)."""
    if not gene_set:
        return annotated_variants, acmg_results
    out_ann: List[Dict[str, Any]] = []
    out_acmg: List[Dict[str, Any]] = []
    for i, a in enumerate(annotated_variants):
        g = (a.get("gene") or "").strip().upper()
        if g not in gene_set:
            continue
        out_ann.append(a)
        if i < len(acmg_results):
            out_acmg.append(acmg_results[i])
        else:
            out_acmg.append({})
    return out_ann, out_acmg


def _realpaths_equal(a: str, b: str) -> bool:
    try:
        return os.path.normpath(os.path.realpath(a)) == os.path.normpath(os.path.realpath(b))
    except OSError:
        return os.path.normpath(os.path.abspath(a)) == os.path.normpath(os.path.abspath(b))


def carrier_result_json_path(job: Job) -> Optional[str]:
    """
    Disk path to ``result.json`` for portal + PDF (must match ``generate_report_json``).

    When ``CARRIER_SCREENING_REPORT_OUTPUT_ROOT`` is set, report generation reads
    ``<root>/output/<work>/<sample>/result.json``, while ``job.output_dir`` may still
    point at the work tree under ``carrier_screening_work_root``. Prefer the same
    directory as ``carrier_report_output_dir`` so dark-genes PATCH and Generate Report
    see one file; fall back to ``job.output_dir`` if only that path exists.

    When ``prior_pipeline_reuse`` left ``job.output_dir`` on the prior order's tree,
    do **not** fall back to that path — it would serve the prior sample's ``result.json``
    as this order's.
    """
    primary = os.path.join(carrier_report_output_dir(job), "result.json")
    if os.path.isfile(primary):
        return primary
    prior_out = (job.params or {}).get("_prior_reuse_output_dir")
    if (
        job.output_dir
        and prior_out
        and (job.params or {}).get("_prior_reuse")
        and _realpaths_equal(str(job.output_dir).strip(), str(prior_out).strip())
    ):
        return None
    if job.output_dir:
        alt = os.path.join(job.output_dir, "result.json")
        if os.path.isfile(alt):
            return alt
    return None


def _atomic_write_text(path: str, text: str, encoding: str = "utf-8") -> None:
    """Write *text* to *path* atomically via temp file + os.replace."""
    import tempfile
    d = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding=encoding) as f:
            f.write(text)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def write_carrier_result_json_sync(job: Job, data: Dict[str, Any]) -> List[str]:
    """
    Write ``result.json`` under ``carrier_report_output_dir`` and mirror to
    ``job.output_dir`` when it differs, so both locations stay aligned.
    Returns paths written.
    """
    primary = os.path.join(carrier_report_output_dir(job), "result.json")
    text = json.dumps(data, ensure_ascii=False, indent=2, default=str)
    os.makedirs(os.path.dirname(primary), exist_ok=True)
    written: List[str] = []
    _atomic_write_text(primary, text)
    written.append(primary)
    if job.output_dir:
        alt = os.path.join(job.output_dir, "result.json")
        if os.path.normpath(alt) != os.path.normpath(primary):
            prior_out = (job.params or {}).get("_prior_reuse_output_dir")
            if (
                (job.params or {}).get("_prior_reuse")
                and prior_out
                and _realpaths_equal(str(job.output_dir).strip(), str(prior_out).strip())
            ):
                pass  # do not mirror into prior order's output tree
            else:
                os.makedirs(os.path.dirname(alt), exist_ok=True)
                _atomic_write_text(alt, text)
                written.append(alt)
    return written


def _project_root_dir() -> str:
    """Repo root (parent of the `app` package)."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))


def _resolve_path_from_settings(p: Optional[str]) -> str:
    """
    Turn env/config paths into absolute paths. Relative values are resolved from the
    repo root (not the process cwd), so REPORT_TEMPLATE_DIR=data/report_templates works
    no matter where uvicorn was started from.
    """
    p = (p or "").strip()
    if not p:
        return ""
    if os.path.isabs(p):
        return os.path.normpath(p)
    return os.path.normpath(os.path.join(_project_root_dir(), p))


def _bundled_carrier_report_template_dir() -> Optional[str]:
    """Shipped templates: <repo>/data/report_templates (Jinja PDF HTML)."""
    d = os.path.join(_project_root_dir(), "data", "report_templates")
    return d if os.path.isdir(d) else None


def _usable_carrier_jinja_dir(path: str) -> bool:
    """Must contain carrier_EN.html or Jinja will fail and we fall back to ugly built-in HTML."""
    if not path or not os.path.isdir(path):
        return False
    return os.path.isfile(os.path.join(path, "carrier_EN.html"))


def _candidate_carrier_template_dirs() -> List[str]:
    """Ordered candidates; first usable (see _usable_carrier_jinja_dir) wins."""
    candidates: List[str] = []
    t = _resolve_path_from_settings(
        getattr(settings, "carrier_screening_report_template_dir", None)
    )
    if t:
        candidates.append(t)
    root = (settings.carrier_screening_report_output_root or "").strip()
    if root:
        root_abs = _resolve_path_from_settings(root)
        # Optional overlay next to pipeline output (prefer new name; legacy folder still supported).
        candidates.append(os.path.join(root_abs, "report_templates"))
        candidates.append(os.path.join(root_abs, "carrier_report"))
    bundled = _bundled_carrier_report_template_dir()
    if bundled:
        candidates.append(bundled)
    legacy = _resolve_path_from_settings(settings.report_template_dir)
    if legacy:
        candidates.append(legacy)
    candidates.append(
        os.path.join(settings.carrier_screening_pipeline_dir, "data", "templates")
    )
    return candidates


def resolve_carrier_pdf_template_dir() -> Optional[str]:
    """
    Directory of Jinja2 templates for WeasyPrint PDFs.

    Only directories that contain carrier_EN.html are used. This avoids picking an
    empty overlay directory before packaged data/report_templates (which would make
    Jinja fail and fall back to built-in HTML).

    Priority:
      1. carrier_screening_report_template_dir
      2. <CARRIER_SCREENING_REPORT_OUTPUT_ROOT>/report_templates
      3. <CARRIER_SCREENING_REPORT_OUTPUT_ROOT>/carrier_report (legacy overlay name)
      4. Packaged data/report_templates
      5. report_template_dir (legacy)
      6. pipeline .../data/templates
    """
    seen: set = set()
    for cand in _candidate_carrier_template_dirs():
        if not cand:
            continue
        norm = os.path.normpath(cand)
        if norm in seen:
            continue
        seen.add(norm)
        if _usable_carrier_jinja_dir(norm):
            logger.info("Carrier PDF: using Jinja template directory %s", norm)
            return norm
        if os.path.isdir(norm):
            logger.warning(
                "Carrier PDF: skipping template directory %r (missing carrier_EN.html)",
                norm,
            )
    logger.error(
        "Carrier PDF: no usable Jinja template directory (need carrier_EN.html). "
        "Reports will use built-in HTML. Add templates under data/report_templates or set "
        "CARRIER_SCREENING_REPORT_TEMPLATE_DIR."
    )
    return None


def _job_stop_requested(order_id: str) -> bool:
    """Portal Stop during PROCESSING: 파이프라인 PID는 없지만 취소 플래그만 설정된 경우."""
    from ...runner import get_runner
    return get_runner().stop_requested(order_id)


def _carrier_requires_interpretation_panel(params: Optional[Dict[str, Any]]) -> bool:
    """
    True when the order must include ``wes_panel_id`` (interpretation/report scope).

    Extended services (portal: ``test_category == other`` + exome/WGS/etc.) hide the
    WES panel selector — those orders skip this requirement.
    """
    p = params or {}
    if not isinstance(p, dict):
        return True
    c = p.get("carrier") or {}
    if not isinstance(c, dict):
        return True
    tc = (c.get("test_category") or "").strip()
    if tc == "other":
        return False
    return True


def _resolve_carrier_wes_panel_id(params: Optional[Dict[str, Any]]) -> str:
    """``wes_panel_id`` at top level or under ``params.carrier`` (portal sends both)."""
    p = params or {}
    if not isinstance(p, dict):
        return ""
    w = (p.get("wes_panel_id") or "").strip()
    if w:
        return w
    c = p.get("carrier") or {}
    if isinstance(c, dict):
        return (c.get("wes_panel_id") or "").strip()
    return ""


class CarrierScreeningPlugin(ServicePlugin):
    """Carrier Screening 서비스 플러그인"""

    # ─── ServicePlugin 메타데이터 ──────────────────────────

    @property
    def service_code(self) -> str:
        return "carrier_screening"

    @property
    def display_name(self) -> str:
        return "Carrier Screening"

    def get_progress_stages(self) -> Dict[str, int]:
        return {
            "ALIGN_BWA": 20,
            "SORT_INDEX": 30,
            "MARK_DUPLICATES": 35,
            "RECALIBRATE": 40,
            "CALL_VARIANTS": 55,
            "GENOTYPE_GVCF": 60,
            "FILTER_VARIANTS": 65,
            "CALL_SV": 70,
            "CALL_CNV": 75,
            "COVERAGE_ANALYSIS": 78,
            "IGV_SNAPSHOT": 80,
            "SUMMARY": 82,
        }

    def validate_params(self, params: Dict[str, Any], strict: bool = True) -> Tuple[bool, str]:
        """서비스별 파라미터 유효성 검사.

        strict=False (Save 단계): wes_panel_id 는 나중에 채울 수 있으므로 건너뜀.
        strict=True  (Submit/Run 단계): 모든 필수 항목 검사.
        """
        carrier = (params or {}).get("carrier") or {}
        if carrier.get("reuse_prior_pipeline_outputs"):
            pid = (carrier.get("prior_order_id") or "").strip()
            if not pid:
                return (
                    False,
                    "prior_order_id is required when reuse_prior_pipeline_outputs is true",
                )
        # Standard carrier screening: full-exome (capture) annotation + order-time interpretation panel.
        # Extended clinical programs (test_category "other" / Exome, WGS, …) do not use the WES panel UI.
        # wes_panel_id 는 Submit/Run 단계(strict=True)에서만 필수 — Save(초안) 단계에서는 선택.
        if strict and _carrier_requires_interpretation_panel(params):
            wid = _resolve_carrier_wes_panel_id(params)
            if not wid:
                return (
                    False,
                    "wes_panel_id is required: choose a Primary (interpretation) panel. "
                    "The pipeline annotates the full exome/capture; the panel narrows the clinical report.",
                )
            from ..wes_panels import get_panel_by_id, resolve_panel_interpretation_genes

            panel = get_panel_by_id(wid)
            if not panel:
                return (
                    False,
                    f"Unknown wes_panel_id={wid!r} — add it to the WES panel catalog or fix the id.",
                )
            genes = resolve_panel_interpretation_genes(panel)
            if not genes:
                return (
                    False,
                    f"Panel {wid!r} resolves to zero interpretation genes "
                    "(needs interpretation_genes, interpretation_genes_file, or "
                    "interpretation_genes_from_disease_bed + disease_bed).",
                )
        return True, ""

    # ─── 디렉토리 구조 ─────────────────────────────────────

    def _get_dirs(self, job: Job) -> Dict[str, str]:
        """
        carrier-screening 디렉토리 구조:
            <work_root>/fastq/<work_segment>/<leaf>/  (work_root = layout_base or CARRIER_SCREENING_ARTIFACT_BASE)

        Use the same ``work_segment`` as ``run_analysis.sh -w`` (``carrier_run_analysis_work_arg``) —
        often from ``.../fastq/<segment>/...`` in R1 — not ``job.work_dir`` when the portal uses a
        different date segment (e.g. work_dir=260403 but FASTQ under fastq/2601/...).
        ``leaf`` is ``carrier_sequencing_folder`` (order_id / explicit sequencing folder).
        """
        work_root = settings.carrier_screening_work_root
        wk = carrier_run_analysis_work_arg(job)
        leaf = carrier_sequencing_folder(job)
        return {
            "fastq": job.fastq_dir
            or os.path.join(work_root, "fastq", wk, leaf),
            "analysis": job.analysis_dir
            or os.path.join(work_root, "analysis", wk, leaf),
            "output": job.output_dir
            or os.path.join(work_root, "output", wk, leaf),
            "log": job.log_dir
            or os.path.join(work_root, "log", wk, leaf),
        }

    def _qc_extra_search_dirs(self, job: Job, analysis_dir: str) -> List[str]:
        """
        artifact 전용 analysis_dir(예: gx-exome-work) 에 QC 파일이 없을 때,
        파이프라인 원본 디렉터리를 추가 검색합니다.

        Nextflow 출력이 job.sample_name 이 아닌 FASTQ 폴더명(sequencing_folder) 아래에만
        있을 수 있으므로, layout_base/analysis/{work}/ 샘플·시퀀싱 폴더 둘 다 후보에 넣습니다.
        """
        wk = carrier_run_analysis_work_arg(job)
        base = settings.carrier_screening_layout_base
        seq = carrier_sequencing_folder(job)
        smp = str(job.sample_name).strip()
        candidates: List[str] = []
        if smp:
            candidates.append(os.path.join(base, "analysis", wk, smp))
        if seq and seq != smp:
            candidates.append(os.path.join(base, "analysis", wk, seq))

        try:
            ad_abs = os.path.abspath(analysis_dir)
        except OSError:
            ad_abs = ""
        seen_real: set = set()
        out: List[str] = []
        for cand in candidates:
            if not cand or not os.path.isdir(cand):
                continue
            try:
                if ad_abs and os.path.abspath(cand) == ad_abs:
                    continue
            except OSError:
                pass
            real = os.path.realpath(cand)
            if real in seen_real:
                continue
            seen_real.add(real)
            out.append(cand)
        return out

    def dark_genes_search_roots(self, job: Job) -> List[str]:
        """
        Same directories as ``generate_result_json`` uses for dark-gene text + visual scans
        (analysis, output, layout extras). Used to merge IGV/repeat paths across roots.
        """
        dirs = self._get_dirs(job)
        analysis_dir = dirs["analysis"]
        output_dir = dirs["output"]
        ad = (analysis_dir or output_dir or "").strip()
        od = (output_dir or "").strip()
        search_roots: List[str] = []
        seen_real: set = set()

        def _add_dark_root(path: Optional[str]) -> None:
            if not path or not str(path).strip():
                return
            p = str(path).strip()
            if not os.path.isdir(p):
                return
            try:
                r = os.path.realpath(p)
            except OSError:
                return
            if r in seen_real:
                return
            seen_real.add(r)
            search_roots.append(os.path.abspath(p))

        try:
            _add_dark_root(ad)
            _add_dark_root(od)
            for extra in self._qc_extra_search_dirs(job, analysis_dir):
                _add_dark_root(extra)
        except OSError:
            pass
        return search_roots

    def _qc_more_roots_from_outputs(self, job: Job) -> List[str]:
        """
        alignment QC(flagstat, MultiQC, Picard)가 analysis 트리가 아니라
        output / layout output / log 에만 있을 때 검색 루트 확장.
        """
        dirs = self._get_dirs(job)
        seen: set = set()
        out: List[str] = []

        def add(p: Optional[str]) -> None:
            if not p:
                return
            p = str(p).strip()
            if not p or not os.path.isdir(p):
                return
            r = os.path.realpath(p)
            if r in seen:
                return
            seen.add(r)
            out.append(p)

        for p in self._vcf_completion_search_roots(job, dirs):
            add(p)
        add(dirs.get("log") or "")
        return out

    def _carrier_bed_directories(self) -> List[str]:
        """파이프라인 패키지 · run_analysis -d 프로젝트 · FASTQ layout 의 data/bed 후보."""
        dirs: List[str] = []
        pd = settings.carrier_screening_pipeline_dir
        dirs.append(os.path.join(pd, "data", "bed"))
        sd = (settings.carrier_screening_script_data_dir or "").strip()
        if sd:
            dirs.append(os.path.join(os.path.abspath(sd), "data", "bed"))
        dirs.append(os.path.join(settings.carrier_screening_layout_base, "data", "bed"))
        seen = set()
        out: List[str] = []
        for d in dirs:
            if d in seen:
                continue
            seen.add(d)
            out.append(d)
        return out

    def _find_default_bed_in_dirs(self, bed_dirs: List[str], prefix: str) -> Optional[str]:
        for bed_dir in bed_dirs:
            if not os.path.isdir(bed_dir):
                continue
            candidates = glob.glob(os.path.join(bed_dir, f"{prefix}*.bed"))
            if candidates:
                return sorted(candidates)[0]
        return None

    def _resolve_disease_gene_json_path(self) -> str:
        candidates: List[str] = []
        if settings.disease_gene_json:
            candidates.append(settings.disease_gene_json)
        pd = settings.carrier_screening_pipeline_dir
        candidates.append(os.path.join(pd, "data", "db", "disease_gene_mapping.json"))
        sd = (settings.carrier_screening_script_data_dir or "").strip()
        if sd:
            candidates.append(
                os.path.join(os.path.abspath(sd), "data", "db", "disease_gene_mapping.json")
            )
        candidates.append(
            os.path.join(settings.carrier_screening_layout_base, "data", "db", "disease_gene_mapping.json")
        )
        seen = set()
        for c in candidates:
            if not c or c in seen:
                continue
            seen.add(c)
            if os.path.isfile(c):
                logger.info("[carrier_screening] Using disease_gene_mapping: %s", c)
                return c
        return settings.disease_gene_json or ""

    def _resolve_carrier_bed(
        self,
        job: Job,
        param_key: str,
        prefix: str,
        default_setting_path: Optional[str],
    ) -> Optional[str]:
        p = job.params.get(param_key)
        if p and os.path.isfile(p):
            return p
        if default_setting_path and os.path.isfile(default_setting_path):
            return default_setting_path
        return self._find_default_bed_in_dirs(self._carrier_bed_directories(), prefix)

    # ─── Step 1: 입력 준비 ─────────────────────────────────

    async def prepare_inputs(self, job: Job) -> bool:
        """FASTQ 파일 확인 및 디렉토리 구조 생성"""
        from .prior_reuse import (
            carrier_reuse_prior_pipeline_requested,
            prepare_carrier_prior_pipeline_reuse,
        )

        dirs = self._get_dirs(job)

        try:
            for d in dirs.values():
                os.makedirs(d, exist_ok=True)
        except PermissionError as e:
            logger.error(
                "[carrier_screening] Permission denied creating %s — fix host ownership, set "
                ".env.compose HOST_UID/HOST_GID to `id -u`/`id -g`, or set CARRIER_SCREENING_ARTIFACT_BASE "
                "to a writable dir under /data (e.g. /data/gx-exome-work) so fastq/analysis/output "
                "staging is not under a read-only carrier layout tree.",
                e.filename or dirs,
            )
            raise

        # Job 경로 업데이트
        job.fastq_dir = dirs["fastq"]
        job.analysis_dir = dirs["analysis"]
        job.output_dir = dirs["output"]
        job.log_dir = dirs["log"]

        # Same sequencing run, new order / panel: reuse prior pipeline outputs (no Nextflow)
        if carrier_reuse_prior_pipeline_requested(job):
            prepare_carrier_prior_pipeline_reuse(job)
            logger.info(
                "[carrier_screening] prior pipeline reuse — skipping FASTQ download/check for %s",
                job.order_id,
            )
            return True

        # Portal / saved order: explicit R1/R2 paths take precedence over directory discovery
        # (avoids wrong pairing when multiple sequencing pairs sit in one folder).
        r1e = (job.fastq_r1_path or "").strip()
        r2e = (job.fastq_r2_path or "").strip()
        if r1e and r2e and os.path.isfile(r1e) and os.path.isfile(r2e):
            perm_err = self._fastq_readable_by_process(r1e, r2e)
            if perm_err:
                logger.error("[carrier_screening] %s", perm_err)
                raise RuntimeError(perm_err)
            job.fastq_r1_path = r1e
            job.fastq_r2_path = r2e
            job.fastq_dir = os.path.dirname(os.path.abspath(r1e))
            logger.info(
                "[carrier_screening] Using explicit FASTQ paths from order: R1=%s",
                r1e,
            )
            return True

        # FASTQ 파일 확인
        r1_path, r2_path = self._find_fastq_files(dirs["fastq"], job.sample_name)

        if r1_path and r2_path:
            perm_err = self._fastq_readable_by_process(r1_path, r2_path)
            if perm_err:
                logger.error("[carrier_screening] %s", perm_err)
                raise RuntimeError(perm_err)
            logger.info(f"Found existing FASTQ files: R1={r1_path}, R2={r2_path}")
            job.fastq_r1_path = r1_path
            job.fastq_r2_path = r2_path
            return True

        # URL에서 다운로드
        if job.fastq_r1_url and job.fastq_r2_url:
            logger.info("Downloading FASTQ files from URLs...")
            r1_dest = os.path.join(dirs["fastq"], f"{job.sample_name}_R1.fastq.gz")
            r2_dest = os.path.join(dirs["fastq"], f"{job.sample_name}_R2.fastq.gz")

            ok1 = await self._download_file(job.fastq_r1_url, r1_dest)
            ok2 = await self._download_file(job.fastq_r2_url, r2_dest)

            if ok1 and ok2:
                perm_err = self._fastq_readable_by_process(r1_dest, r2_dest)
                if perm_err:
                    logger.error("[carrier_screening] %s", perm_err)
                    raise RuntimeError(perm_err)
                job.fastq_r1_path = r1_dest
                job.fastq_r2_path = r2_dest
                return True
            else:
                msg = "Failed to download FASTQ files from URLs (wget returned non-zero)"
                logger.error("[carrier_screening] %s", msg)
                raise RuntimeError(msg)

        # 지정된 로컬 경로 확인 (exists but not isfile — e.g. broken symlink)
        if job.fastq_r1_path and job.fastq_r2_path:
            if os.path.isfile(job.fastq_r1_path) and os.path.isfile(job.fastq_r2_path):
                perm_err = self._fastq_readable_by_process(
                    job.fastq_r1_path, job.fastq_r2_path
                )
                if perm_err:
                    logger.error("[carrier_screening] %s", perm_err)
                    raise RuntimeError(perm_err)
                job.fastq_dir = os.path.dirname(os.path.abspath(job.fastq_r1_path))
                return True
            r1_ok = os.path.isfile((job.fastq_r1_path or "").strip())
            r2_ok = os.path.isfile((job.fastq_r2_path or "").strip())
            msg = (
                f"FASTQ paths set on the order but not both readable files "
                f"(R1 ok={r1_ok} R2 ok={r2_ok}): "
                f"R1={job.fastq_r1_path!r} R2={job.fastq_r2_path!r}"
            )
            logger.error("[carrier_screening] %s", msg)
            raise RuntimeError(msg)

        msg = (
            f"No FASTQ files found for sample_name={job.sample_name!r} under {dirs['fastq']!r}. "
            f"Set R1/R2 paths on the order, place *_R1_* / *_R2_* (or _1._/_2.) .fastq.gz in that folder, "
            f"or enable prior pipeline reuse with a completed carrier/exome order."
        )
        logger.error("[carrier_screening] %s", msg)
        raise RuntimeError(msg)

    @staticmethod
    def _fastq_readable_by_process(r1: str, r2: str) -> Optional[str]:
        """
        Return an error message if FASTQs exist but cannot be read by this process.

        Typical host issue: files are mode 600/700 and owned by another user — bwa in the
        Nextflow task then fails with 'fail to open file' even when paths are correct.
        """
        import stat as stat_mod

        for label, p in (("R1", r1), ("R2", r2)):
            if not p or not str(p).strip():
                return f"{label} FASTQ path is empty"
            if not os.path.isfile(p):
                return f"{label} FASTQ not found: {p}"
            if not os.access(p, os.R_OK):
                try:
                    st = os.stat(p)
                    mode = stat_mod.filemode(st.st_mode)
                except OSError:
                    mode = "?"
                return (
                    f"{label} FASTQ not readable by the daemon user ({mode}): {p}. "
                    "On the host: chmod a+r (or chown to HOST_UID from .env.compose), "
                    "or chmod g+rX and run the daemon as a user in the file's group."
                )
        return None

    def _find_fastq_files(self, fastq_dir: str, sample_name: str) -> Tuple[Optional[str], Optional[str]]:
        """
        FASTQ 디렉토리에서 R1/R2 파일을 찾습니다.

        파일명 패턴:
            - {sample}_R1_*.fastq.gz / {sample}_R2_*.fastq.gz
            - {sample}_1.fastq.gz / {sample}_2.fastq.gz
            - {sample}_R1.fq.gz / {sample}_R2.fq.gz
        """
        if not os.path.isdir(fastq_dir):
            return None, None

        all_files = os.listdir(fastq_dir)
        fastq_files = [
            f for f in all_files
            if f.endswith((".fastq.gz", ".fq.gz"))
        ]

        if not fastq_files:
            return None, None

        r1 = None
        r2 = None

        # 패턴 1: _R1_ / _R2_
        for f in fastq_files:
            if re.search(r"_R1[_.]", f, re.IGNORECASE):
                r1 = os.path.join(fastq_dir, f)
            elif re.search(r"_R2[_.]", f, re.IGNORECASE):
                r2 = os.path.join(fastq_dir, f)

        # 패턴 2: _1. / _2.
        if not r1 or not r2:
            for f in fastq_files:
                if re.search(r"_1\.(fastq|fq)\.gz$", f, re.IGNORECASE):
                    r1 = r1 or os.path.join(fastq_dir, f)
                elif re.search(r"_2\.(fastq|fq)\.gz$", f, re.IGNORECASE):
                    r2 = r2 or os.path.join(fastq_dir, f)

        return r1, r2

    async def _download_file(self, url: str, dest: str) -> bool:
        """URL에서 파일을 다운로드합니다."""
        try:
            proc = await asyncio.create_subprocess_exec(
                "wget", "-q", "-O", dest, url,
                stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE,
            )
            _, stderr = await proc.communicate()
            if proc.returncode != 0:
                logger.error(f"Download failed: {stderr.decode()}")
                return False
            return os.path.exists(dest) and os.path.getsize(dest) > 0
        except Exception as e:
            logger.error(f"Download error: {e}")
            return False

    # ─── Step 2: 파이프라인 실행 명령 ──────────────────────

    def _effective_run_analysis_script(self) -> Optional[str]:
        """
        1) CARRIER_SCREENING_RUN_SCRIPT (절대 경로)
        2) 프로젝트 루트(layout base, 예: /data/gx-exome) 아래 src/run_analysis.sh
           — compose 가 CARRIER_SCREENING_HOST 로 레포 전체를 붙인 경우
        3) CARRIER_SCREENING_PIPELINE_DIR 아래 src/, 루트, scripts/
        없으면 None → Nextflow 직접.

        CARRIER_SCREENING_FORCE_DIRECT_NEXTFLOW=1 skips 1–3 so DARK_GENE_PIPELINE_DIR can be used.
        """
        if settings.carrier_screening_force_direct_nextflow:
            logger.info(
                "[carrier_screening] CARRIER_SCREENING_FORCE_DIRECT_NEXTFLOW set — "
                "skipping run_analysis.sh, using direct Nextflow"
            )
            return None
        explicit = (settings.carrier_screening_run_script or "").strip()
        layout_base = settings.carrier_screening_layout_base
        pd = settings.carrier_screening_pipeline_dir
        to_try: List[str] = []
        if explicit:
            to_try.append(explicit)
        for rel in ("src/run_analysis.sh", "run_analysis.sh", "scripts/run_analysis.sh"):
            p = os.path.join(layout_base, rel)
            if p not in to_try:
                to_try.append(p)
        for rel in ("src/run_analysis.sh", "run_analysis.sh", "scripts/run_analysis.sh"):
            p = os.path.join(pd, rel)
            if p not in to_try:
                to_try.append(p)
        for path in to_try:
            if path and os.path.isfile(path):
                if not explicit or path != explicit:
                    logger.info("[carrier_screening] Using run_analysis.sh at %s", path)
                return path
        if explicit:
            logger.warning(
                "[carrier_screening] CARRIER_SCREENING_RUN_SCRIPT not found: %s; "
                "also checked %s and %s — using Nextflow direct",
                explicit,
                layout_base,
                pd,
            )
        else:
            logger.warning(
                "[carrier_screening] No run_analysis.sh — set CARRIER_SCREENING_RUN_SCRIPT or mount "
                "repo at layout base (e.g. %s/src/run_analysis.sh) or under %s — using Nextflow direct",
                layout_base,
                pd,
            )
        return None

    def _run_analysis_data_dir(self) -> str:
        """
        run_analysis.sh --data-dir: must be a **host** path for nested `docker run -v`.

        Prefer CARRIER_SCREENING_SCRIPT_DATA_DIR. If unset and layout is the container mount
        /data/gx-exome, map via CARRIER_SCREENING_HOST (same as compose) so Docker binds real dirs.
        """
        d = (settings.carrier_screening_script_data_dir or "").strip()
        if d:
            return d
        host = os.environ.get("CARRIER_SCREENING_HOST", "").strip()
        layout = settings.carrier_screening_layout_base
        if host and layout.startswith(_GX_EXOME_LAYOUT_PREFIX):
            suffix = layout[len(_GX_EXOME_LAYOUT_PREFIX) :].lstrip("/")
            if suffix:
                return os.path.join(host.rstrip("/"), suffix)
            return host
        return layout

    def _run_analysis_ref_dir(self) -> str:
        r = (settings.carrier_screening_script_ref_dir or "").strip()
        if r:
            return r
        return os.path.join(self._run_analysis_data_dir(), "data", "refs")

    def _ensure_fastq_symlinks(self, job: Job, expected_fastq_dir: str) -> None:
        """
        fastq_r1_path / fastq_r2_path 가 expected_fastq_dir 밖에 있을 때,
        expected_fastq_dir 안에 상대 심볼릭 링크를 생성한다.

        목적:
          - --sample = order_id 로 고정하되, FASTQ 는 어느 디렉토리에서든 선택 가능.
          - 상대 링크이므로 docker DinD 에서 ${DATA_DIR}/fastq 를 마운트해도 정상 작동.
          - 읽기 전용 FASTQ 를 여러 오더가 동시에 사용해도 안전 (read-only, 상태 없음).
        """
        for attr in ("fastq_r1_path", "fastq_r2_path"):
            fq = (getattr(job, attr, None) or "").strip()
            if not fq:
                continue
            fq_abs = os.path.abspath(fq)
            if not os.path.isfile(fq_abs):
                continue
            fq_dir = os.path.realpath(os.path.dirname(fq_abs))
            expected_real = os.path.realpath(expected_fastq_dir) if os.path.isdir(expected_fastq_dir) else expected_fastq_dir
            if fq_dir == expected_real:
                continue  # 이미 올바른 위치

            os.makedirs(expected_fastq_dir, exist_ok=True)
            fname = os.path.basename(fq_abs)
            link_path = os.path.join(expected_fastq_dir, fname)
            rel_target = os.path.relpath(fq_abs, expected_fastq_dir)

            if os.path.islink(link_path):
                existing_target = os.readlink(link_path)
                if existing_target == rel_target:
                    continue  # 올바른 링크 — 건너뜀
                logger.info(
                    "[carrier_screening] Replacing stale FASTQ symlink: %s (was %s → now %s)",
                    link_path, existing_target, rel_target,
                )
                os.unlink(link_path)
            elif os.path.exists(link_path):
                continue  # 실제 파일 — 건드리지 않음

            try:
                os.symlink(rel_target, link_path)
                logger.info(
                    "[carrier_screening] FASTQ symlink created: %s → %s (order=%s)",
                    link_path, rel_target, job.order_id,
                )
            except OSError as e:
                logger.warning(
                    "[carrier_screening] FASTQ symlink failed: %s → %s: %s — FASTQ must be copied manually.",
                    link_path, fq_abs, e,
                )

    def _resolve_capture_panel_bed(self, capture_panel_id: str) -> Optional[str]:
        """
        capture_panel_id → run_analysis.sh 에 넘길 --backbone-bed 절대 경로.

        탐색 순서:
          1. CARRIER_CAPTURE_PANEL_BED_DIR/{capture_panel_id}/targets.bed  (명시적 설정)
          2. {CARRIER_SCREENING_SCRIPT_DATA_DIR}/data/bed/{capture_panel_id}/targets.bed  (기본)

        반환값은 HOST 경로 (run_analysis.sh 내부 docker run에서 직접 참조).
        존재하지 않으면 None → --panel 폴백.
        """
        data_dir = self._run_analysis_data_dir()
        candidates: List[str] = []

        explicit_base = (settings.carrier_capture_panel_bed_dir or "").strip()
        if explicit_base:
            candidates.append(os.path.join(explicit_base, capture_panel_id, "targets.bed"))

        # 기본: DATA_DIR/data/bed/{panel}/targets.bed  (gx-exome 표준 레이아웃)
        candidates.append(os.path.join(data_dir, "data", "bed", capture_panel_id, "targets.bed"))

        for path in candidates:
            if os.path.isfile(path):
                logger.info(
                    "[carrier_screening] capture panel BED resolved: panel=%s path=%s",
                    capture_panel_id, path,
                )
                return path

        logger.warning(
            "[carrier_screening] capture panel BED not found for panel=%s (tried: %s) — falling back to --panel",
            capture_panel_id, candidates,
        )
        return None

    def _shell_command_run_analysis(self, job: Job, script: str) -> str:
        """
        gx-exome run_analysis.sh 호출 (long-form 인자).
        예: bash run_analysis.sh --work-dir 2601 --sample SAMPLE --data-dir /home/ken/gx-exome
                --backbone-bed /home/ken/gx-exome/data/bed/twist-exome2/targets.bed
                --aligner bwa-mem2 --variant-caller deepvariant --no-skip-vep --skip-cnv

        --backbone-bed 우선: service-daemon 이 패널→BED 경로를 직접 관리.
        BED 파일이 없을 경우 --panel capture_panel_id 로 폴백.
        --sample = order_id (항상 고정).
          FASTQ 파일이 다른 디렉토리에 있을 경우 {data_dir}/fastq/{work}/{order_id}/ 아래에
          상대 심볼릭 링크를 자동 생성하여 run_analysis.sh 에 올바른 --sample 을 전달한다.
          output/analysis/log 는 항상 order_id 기준으로 생성된다.
        """
        data_dir = self._run_analysis_data_dir()
        sample_folder = carrier_sequencing_folder(job)   # == order_id (명시 재정의 없는 한)
        work_arg = carrier_run_analysis_work_arg(job)
        carrier_params = (job.params or {}).get("carrier") or {}
        capture_panel = (carrier_params.get("capture_panel_id") or "").strip() or "twist-exome2"

        # FASTQ가 예상 위치({data_dir}/fastq/{work}/{order_id}/)에 없으면 심볼릭 링크 생성.
        # 같은 fastq/ 트리 내의 상대 링크이므로 docker DinD 마운트에서도 정상 작동.
        fastq_dir = os.path.join(data_dir, "fastq", work_arg, sample_folder)
        self._ensure_fastq_symlinks(job, fastq_dir)

        if not os.path.isdir(fastq_dir):
            raise FileNotFoundError(
                f"[carrier_screening] FASTQ directory not found: {fastq_dir}\n"
                f"  order_id='{job.order_id}', --sample='{sample_folder}'\n"
                f"  FASTQ를 {fastq_dir}/ 에 위치시키거나 fastq_r1_path/r2_path 를 지정하세요."
            )

        parts = [
            _bash_executable(),
            script,
            "--work-dir", work_arg,
            "--sample", sample_folder,
            "--data-dir", data_dir,
        ]

        # service-daemon 이 BED 경로를 알면 --backbone-bed 로 직접 전달 (패널 관리 주체 명확화)
        bed_path = self._resolve_capture_panel_bed(capture_panel)
        if bed_path:
            parts += ["--backbone-bed", bed_path]
        else:
            parts += ["--panel", capture_panel]

        parts += [
            "--aligner", "bwa-mem2",
            "--variant-caller", "deepvariant",
            "--no-skip-vep",
            "--skip-cnv",
        ]

        is_fresh = bool((job.params or {}).get("_pipeline_fresh"))
        logger.info(
            "[carrier_screening] _shell_command_run_analysis: order=%s sample=%s panel=%s bed=%s fresh=%s",
            job.order_id, sample_folder, capture_panel, bed_path or "(--panel fallback)", is_fresh,
        )
        if is_fresh:
            parts.append("--fresh")
        extra = (settings.carrier_screening_script_extra_args or "").strip()
        if extra:
            parts.extend(shlex.split(extra))
        return " ".join(shlex.quote(p) for p in parts)

    def _resolve_nextflow_main_and_config(
        self,
    ) -> Tuple[str, Optional[str]]:
        """
        main.nf / nextflow.config 위치 — 레포마다 bin/, 루트, workflow/ 등 다름.
        """
        pd = settings.carrier_screening_pipeline_dir.rstrip("/")
        explicit = (settings.carrier_screening_main_nf or "").strip()
        if explicit:
            main_nf = os.path.abspath(explicit)
            if os.path.isfile(main_nf):
                mdir = os.path.dirname(main_nf)
                for cand in (
                    os.path.join(mdir, "nextflow.config"),
                    os.path.join(pd, "nextflow.config"),
                    os.path.join(pd, "bin", "nextflow.config"),
                ):
                    if os.path.isfile(cand):
                        return main_nf, cand
                return main_nf, None
            logger.warning(
                "[carrier_screening] CARRIER_SCREENING_MAIN_NF not found: %s — scanning pipeline_dir",
                explicit,
            )

        dg = (settings.dark_gene_pipeline_dir or "").strip()
        if dg:
            dg_main = os.path.join(os.path.abspath(dg), "main.nf")
            if os.path.isfile(dg_main):
                mdir = os.path.dirname(dg_main)
                for cand in (
                    os.path.join(mdir, "nextflow.config"),
                    os.path.join(pd, "nextflow.config"),
                ):
                    if os.path.isfile(cand):
                        return dg_main, cand
                return dg_main, None

        pairs = [
            ("bin/main.nf", "bin/nextflow.config"),
            ("main.nf", "nextflow.config"),
            ("workflow/main.nf", "workflow/nextflow.config"),
            ("workflow/main.nf", "nextflow.config"),
            ("workflows/main.nf", "nextflow.config"),
        ]
        tried: List[str] = []
        for main_rel, cfg_rel in pairs:
            main_path = os.path.join(pd, main_rel)
            tried.append(main_path)
            if not os.path.isfile(main_path):
                continue
            cfg_path = os.path.join(pd, cfg_rel)
            if os.path.isfile(cfg_path):
                return main_path, cfg_path
            side = os.path.join(os.path.dirname(main_path), "nextflow.config")
            if os.path.isfile(side):
                return main_path, side
            root_cfg = os.path.join(pd, "nextflow.config")
            if os.path.isfile(root_cfg):
                return main_path, root_cfg
            return main_path, None

        fallback = os.path.join(pd, "bin", "main.nf")
        logger.error(
            "[carrier_screening] No main.nf found under %s (tried: %s). "
            "Set CARRIER_SCREENING_MAIN_NF to the real path.",
            pd,
            ", ".join(tried),
        )
        return fallback, None

    def _uses_dark_gene_unified(self, main_nf: str) -> bool:
        """True when main.nf is the dark_gene_pipeline single workflow (expects --fastq_dir, ref_fasta, outdir)."""
        dg = (settings.dark_gene_pipeline_dir or "").strip()
        if dg:
            exp = os.path.abspath(os.path.join(dg, "main.nf"))
            try:
                if os.path.abspath(main_nf) == exp:
                    return True
            except OSError:
                pass
        low = (main_nf or "").lower()
        if "dark_gene" in low and low.endswith("main.nf"):
            return True
        return False

    def _dark_gene_unified_params(self, job: Job, main_nf: str) -> Dict[str, Any]:
        """CLI params for dark_gene_pipeline/main.nf (see that repo nextflow.config)."""
        r1 = (job.fastq_r1_path or "").strip()
        if not r1:
            raise ValueError("fastq_r1_path required for dark gene unified pipeline")
        fastq_dir = os.path.dirname(os.path.abspath(r1))
        root = os.path.dirname(os.path.abspath(main_nf))
        out: Dict[str, Any] = {
            "fastq_dir": fastq_dir,
            "outdir": job.analysis_dir,
        }
        if settings.ref_fasta:
            out["ref_fasta"] = settings.ref_fasta
        if settings.ref_fai:
            out["ref_fai"] = settings.ref_fai
        if settings.ref_dict:
            out["ref_dict"] = settings.ref_dict
        if settings.ref_bwa_indices:
            out["ref_bwa_indices"] = settings.ref_bwa_indices
        backbone = self._resolve_carrier_bed(
            job,
            "backbone_bed",
            "backbone",
            settings.carrier_default_backbone_bed,
        )
        if backbone:
            out["backbone_bed"] = backbone
        eh = os.path.join(root, "bed", "variant_catalog_grch38.json")
        if os.path.isfile(eh):
            out["eh_catalog"] = eh
        out["skip_cnv"] = settings.dark_gene_skip_cnv
        for key in ("pon_tar", "gcnv_model", "interval_list", "hba_bed", "cyp21a2_bed", "eh_catalog", "cleanup"):
            v = job.params.get(key)
            if v:
                out[key] = v
        return out

    async def get_pipeline_command(self, job: Job) -> str:
        """src/run_analysis.sh 가 있으면 우선 사용, 없으면 Nextflow 직접 실행."""
        from ..wes_panels import apply_wes_panel_to_job_params

        apply_wes_panel_to_job_params(job)

        if job.params.get("_prior_reuse"):
            logger.info(
                "[carrier_screening] prior pipeline reuse — skipping Nextflow/run_analysis for %s",
                job.order_id,
            )
            return "true"

        script_path = self._effective_run_analysis_script()
        if script_path:
            if not os.path.isfile(script_path):
                raise RuntimeError(
                    f"Carrier run script not found: {script_path!r} "
                    "(mount gx-exome at /data/gx-exome or set CARRIER_SCREENING_RUN_SCRIPT)"
                )
            cmd = self._shell_command_run_analysis(job, script_path)
            logger.info(
                "[carrier_screening] Using run script %s: %s",
                script_path,
                cmd if len(cmd) <= 240 else cmd[:240] + "...",
            )
            return cmd

        logger.warning(
            "[carrier_screening] Pipeline command will use Nextflow (see previous warning if run_analysis.sh missing)"
        )
        pipeline_dir = settings.carrier_screening_pipeline_dir
        main_nf, nf_config = self._resolve_nextflow_main_and_config()
        logger.info(
            "[carrier_screening] Nextflow entry: %s (config: %s)",
            main_nf,
            nf_config or "(none)",
        )

        # Nextflow 실행 명령 구성 (절대 경로로 127 방지)
        cmd_parts = [
            _resolve_nextflow_executable(),
            "run",
            main_nf,
        ]

        prof = (settings.carrier_screening_nextflow_profile or "").strip()
        if prof:
            cmd_parts.extend(["-profile", prof])

        # config 파일
        if nf_config and os.path.isfile(nf_config):
            cmd_parts.extend(["-c", nf_config])

        use_dark = self._uses_dark_gene_unified(main_nf)
        params: Dict[str, Any]
        if use_dark and (job.fastq_r1_path or "").strip():
            try:
                params = self._dark_gene_unified_params(job, main_nf)
                logger.info(
                    "[carrier_screening] Dark gene unified Nextflow (fastq_dir=%s, outdir=%s)",
                    params.get("fastq_dir"),
                    params.get("outdir"),
                )
            except ValueError as e:
                logger.error("[carrier_screening] %s — falling back to legacy NF params", e)
                use_dark = False
        if not use_dark or not (job.fastq_r1_path or "").strip():
            if use_dark and not (job.fastq_r1_path or "").strip():
                logger.warning(
                    "[carrier_screening] Dark gene unified selected but fastq_r1_path missing; "
                    "using legacy carrier params (Nextflow may fail)."
                )
            params = {
                "sample_name": carrier_sequencing_folder(job),
                "fastq_r1": job.fastq_r1_path,
                "fastq_r2": job.fastq_r2_path,
                "outdir": job.analysis_dir,
            }
            if settings.ref_fasta:
                params["ref"] = settings.ref_fasta
            if settings.ref_bwa_indices:
                params["bwa_index"] = settings.ref_bwa_indices
            backbone_bed = self._resolve_carrier_bed(
                job,
                "backbone_bed",
                "backbone",
                settings.carrier_default_backbone_bed,
            )
            if backbone_bed:
                params["backbone_bed"] = backbone_bed
            for key in ("pon_tar", "target_bed", "disease_bed", "cnv_bed"):
                if job.params.get(key):
                    params[key] = job.params[key]

        # 파라미터를 명령줄에 추가 (omit None / empty string; omit False except skip_cnv)
        for key, val in params.items():
            if val is None or val == "":
                continue
            if val is False and key != "skip_cnv":
                continue
            cmd_parts.append(f"--{key}")
            if isinstance(val, bool):
                cmd_parts.append("true" if val else "false")
            else:
                cmd_parts.append(str(val))

        # 작업 디렉토리
        work_dir = os.path.join(job.analysis_dir, "work")
        cmd_parts.extend(["-work-dir", work_dir])

        # fresh 모드: work/ 와 .nextflow 세션 캐시를 삭제하고 -resume 없이 실행
        if (job.params or {}).get("_pipeline_fresh"):
            import shutil
            nf_cache = os.path.join(job.analysis_dir, ".nextflow")
            for d in (work_dir, nf_cache):
                if os.path.isdir(d):
                    logger.info("[carrier_screening] --fresh: removing %s", d)
                    shutil.rmtree(d, ignore_errors=True)
        else:
            cmd_parts.append("-resume")

        # 리포트
        report_path = os.path.join(job.log_dir, "nextflow_report.html")
        cmd_parts.extend(["-with-report", report_path])

        return shlex.join(cmd_parts)

    # ─── Step 3: 완료 확인 ─────────────────────────────────

    def _vcf_completion_search_roots(self, job: Job, dirs: Dict[str, str]) -> List[str]:
        """
        VCF 후보 경로:
        - Job이 가리키는 analysis / output (보통 CARRIER_SCREENING_ARTIFACT_BASE 아래)
        - carrier_screening_layout_base/output/<work>/... (FASTQ 프로젝트 루트)
        - carrier_screening_script_data_dir/output/<work>/... (run_analysis -d 트리; 호스트 경로)

        ARTIFACT_BASE 를 쓰면 Nextflow outdir 은 work 아래인데, run_analysis 는 프로젝트 루트
        아래 output/.../vcf 에만 쓰는 경우가 많다.
        """
        roots: List[str] = []

        def add(path: str) -> None:
            p = (path or "").strip()
            if p and os.path.isdir(p):
                ap = os.path.abspath(p)
                if ap not in roots:
                    roots.append(ap)

        for key in ("_prior_reuse_analysis_dir", "_prior_reuse_output_dir"):
            add((job.params or {}).get(key) or "")

        for key in ("analysis", "output"):
            add(dirs.get(key) or "")

        layout_base = settings.carrier_screening_layout_base
        work_root = os.path.abspath(settings.carrier_screening_work_root)
        script_root_raw = (settings.carrier_screening_script_data_dir or "").strip()
        # Pipeline -w / output tree follows FASTQ .../fastq/<segment>/... which may differ from job.work_dir
        # (e.g. work_dir=260403 for the order but FASTQ under fastq/2601/...).
        wk_run = carrier_run_analysis_work_arg(job)
        wk_job = str(job.work_dir).strip() or "00"
        work_segments: List[str] = []
        for w in (wk_run, wk_job):
            if w and w not in work_segments:
                work_segments.append(w)
        seq = carrier_sequencing_folder(job)
        path_key = (str(job.sample_name).strip() or job.order_id or "").strip()
        leaves = list(dict.fromkeys([x for x in (path_key, seq) if x]))

        seen_bases = set()
        for base in (layout_base, script_root_raw):
            if not (base or "").strip():
                continue
            ba = os.path.abspath(base.strip())
            if ba in seen_bases:
                continue
            seen_bases.add(ba)
            # Per-leaf output|analysis under layout/script roots. Do **not** skip when ba == work_root:
            # job dirs use sample_name only, but run_analysis writes under sequencing_folder when it
            # differs from sample_name — skipping the whole block hid output/<work>/<seq>/vcf/.
            for wseg in work_segments:
                for leaf in leaves:
                    add(os.path.join(ba, "output", wseg, leaf))
                    add(os.path.join(ba, "analysis", wseg, leaf))

        return roots

    def _carrier_vcf_basename_tokens(self, job: Job) -> List[str]:
        """파이프 산출 파일명이 sample_name 전체와 다를 때(예: NA12878만 포함) 대비."""
        seq = carrier_sequencing_folder(job)
        raw: List[Optional[str]] = [job.sample_name, seq, job.order_id]
        if job.params.get("_prior_reuse"):
            raw.extend(
                [
                    job.params.get("_prior_reuse_sample_name"),
                    job.params.get("_prior_reuse_order_id"),
                ]
            )
        tokens: List[str] = []
        for r in raw:
            if not r:
                continue
            s = str(r).strip()
            if not s:
                continue
            tokens.append(s)
            for part in re.split(r"[-_]", s):
                p = part.strip()
                if len(p) >= 4:
                    tokens.append(p)
        seen = set()
        out: List[str] = []
        for t in tokens:
            k = t.lower()
            if k not in seen:
                seen.add(k)
                out.append(t)
        return out

    def _carrier_main_vcf_sort_key(self, path: str, tokens: List[str]) -> Tuple:
        """
        소형 변이용 메인 VCF 우선: 이름에 *annotated* (VEP CSQ 등) > *filtered* (주석 없을 수 있음)
        > *genotyped* > 기타. SV(manta)·repeat 등은 후순위.
        """
        norm = path.replace("\\", "/")
        base = os.path.basename(norm)
        bl = base.lower()
        ll = norm.lower()
        in_vcf = "/vcf/" in ll or ll.rstrip("/").endswith("/vcf")
        in_sv = "/sv/" in ll
        in_rep = "/repeat/" in ll
        has_ann = "annotated" in bl   # VEP CSQ 포함 → 최우선
        has_filt = "filtered" in bl and not has_ann  # annotation 없는 filtered만 후순위
        has_gen = "genotyped" in bl
        has_manta = "manta" in bl
        token_hit = any(t.lower() in bl for t in tokens)
        # annotated(0) > filtered(1) > genotyped(2) > 기타(3)
        ann_rank = 0 if has_ann else (1 if has_filt else (2 if has_gen else 3))
        return (
            0 if in_vcf else 1,
            ann_rank,
            1 if (in_sv or in_rep) else 0,
            1 if has_manta else 0,
            0 if token_hit else 1,
            base,
        )

    def _pick_main_carrier_vcf(self, vcf_files: List[str], job: Job) -> str:
        tokens = self._carrier_vcf_basename_tokens(job)
        return min(vcf_files, key=lambda p: self._carrier_main_vcf_sort_key(p, tokens))

    async def check_completion(self, job: Job) -> bool:
        """파이프라인 완료를 확인합니다 (VCF 파일 존재 여부)."""
        dirs = self._get_dirs(job)
        search_roots = self._vcf_completion_search_roots(job, dirs)

        vcf_files: List[str] = []
        for root in search_roots:
            # depth 3 이내에서만 탐색 (output/<work>/<sample>/vcf/ 구조)
            for depth_pat in ("*.vcf.gz", "*/*.vcf.gz", "*/*/*.vcf.gz", "*/*/*/*.vcf.gz"):
                vcf_files.extend(glob.glob(os.path.join(root, depth_pat)))
            for depth_pat in ("*.vcf", "*/*.vcf", "*/*/*.vcf", "*/*/*/*.vcf"):
                vcf_files.extend(glob.glob(os.path.join(root, depth_pat)))
        vcf_files = list(dict.fromkeys(vcf_files))

        if not vcf_files:
            logger.error(
                "[carrier_screening] No VCF under any search root (artifact analysis/output, "
                "layout_base/output, script_data_dir/output): %s",
                " | ".join(search_roots) if search_roots else "(none)",
            )
            return False

        hint_raw = (job.params or {}).get("_prior_reuse_main_vcf_hint")
        hint = ""
        if isinstance(hint_raw, str) and hint_raw.strip():
            hint = normalize_legacy_carrier_container_path(hint_raw.strip()) or hint_raw.strip()
        main_vcf: str
        if hint and os.path.isfile(hint):
            by_real: Dict[str, str] = {}
            for p in vcf_files:
                try:
                    by_real[os.path.realpath(os.path.abspath(p))] = p
                except OSError:
                    by_real[p] = p
            try:
                hreal = os.path.realpath(os.path.abspath(hint))
            except OSError:
                hreal = hint
            if hreal in by_real:
                main_vcf = by_real[hreal]
            else:
                main_vcf = self._pick_main_carrier_vcf(vcf_files, job)
        else:
            main_vcf = self._pick_main_carrier_vcf(vcf_files, job)
        logger.info(
            "[carrier_screening] Main VCF (ranked among %s candidates): %s",
            len(vcf_files),
            main_vcf,
        )
        job.params["main_vcf"] = normalize_legacy_carrier_container_path(main_vcf) or main_vcf
        return True

    # ─── Step 4: 결과 후처리 (Annotation) ──────────────────

    async def process_results(self, job: Job) -> bool:
        """
        VCF → BED 필터 → Annotation → ACMG 분류 → result.json 생성

        통합 parse_vcf_variants 파이프라인을 사용하여 처리합니다.

        전체 annotation 파이프라인:
            1. VCF FORMAT 문제 필드 정리 (이미 main보다 새 cleaned 가 있으면 생략)
            2. snpEff annotation (선택적, cleaned 이상으로 최신 snpEff 가 있으면 생략)
            3. VariantAnnotator 초기화 (모든 DB 소스 포함)
            4. VariantFilterConfig 구성 (BED, HPO, AF, ClinVar 등)
            5. 통합 parse_vcf_variants 실행
            6. QC 메트릭스 추출
            7. result.json + variants.tsv 생성
        """
        try:
            from .vcf_parser import (
                clean_vcf_remove_formats, run_snpeff,
                load_bed_regions, parse_vcf_variants,
                VariantFilterConfig,
            )
            from .annotator import VariantAnnotator
            from .review import extract_qc_summary, generate_result_json, generate_variants_tsv
            from .vep_parser import is_vep_annotated_vcf, extract_vep_annotations_from_vcf
            from .layout_norm import align_carrier_job_dirs_from_main_vcf

            # ── 0. main VCF 실제 위치와 job.analysis|output|log 정렬 (work_root ≠ 파이프 산출 트리인 경우)
            main_raw = job.params.get("main_vcf")
            main_vcf = (
                normalize_legacy_carrier_container_path(main_raw)
                if isinstance(main_raw, str)
                else None
            )
            if main_vcf and main_raw != main_vcf:
                job.params["main_vcf"] = main_vcf
            if not main_vcf or not os.path.exists(main_vcf):
                logger.error(f"Main VCF not found: {main_vcf}")
                return False
            align_carrier_job_dirs_from_main_vcf(job, main_vcf)

            from ..wes_panels import apply_wes_panel_to_job_params

            apply_wes_panel_to_job_params(job)

            dirs = self._get_dirs(job)
            analysis_dir = dirs["analysis"]
            # Pipeline artifacts (VCF, pipeline_complete) may live under a prior run's tree
            # after main_vcf alignment; portal outputs go to carrier_report_output_dir (this order).
            artifact_output_dir = dirs["output"]
            os.makedirs(artifact_output_dir, exist_ok=True)
            portal_output_dir = carrier_report_output_dir(job)
            os.makedirs(portal_output_dir, exist_ok=True)

            logger.info(f"[process_results] Starting annotation for {job.sample_name}")
            logger.info(f"  Main VCF: {main_vcf}")
            logger.info(
                "[process_results] artifact output (VCF tree)=%s | portal result.json dir=%s",
                artifact_output_dir,
                portal_output_dir,
            )
            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested before cleaning VCF")
                return False

            # ── 1b. VEP annotation 여부 감지 ──
            # Nextflow 파이프라인에서 VEP annotation이 완료된 경우
            # (pipeline_complete.json의 vep_annotation == 'enabled' 또는 VCF 헤더 직접 확인)
            pipeline_complete_path = os.path.join(artifact_output_dir, "pipeline_complete.json")
            prior_od = (job.params or {}).get("_prior_reuse_output_dir")
            if not os.path.exists(pipeline_complete_path) and prior_od:
                alt_pc = os.path.join(str(prior_od).strip(), "pipeline_complete.json")
                if os.path.isfile(alt_pc):
                    pipeline_complete_path = alt_pc
            vep_enabled_from_pipeline = False
            if os.path.exists(pipeline_complete_path):
                try:
                    import json as _json
                    with open(pipeline_complete_path) as _f:
                        _pc = _json.load(_f)
                    vep_enabled_from_pipeline = _pc.get("vep_annotation", "") == "enabled"
                except Exception:
                    pass

            has_csq_header = is_vep_annotated_vcf(main_vcf)
            is_vep_vcf = vep_enabled_from_pipeline or has_csq_header
            logger.info(f"  VEP annotation detected: {is_vep_vcf}")
            logger.info(
                "[process_results] main_vcf fingerprint — path=%s basename=%s "
                "##INFO_CSQ_header=%s pipeline_complete.vep_annotation_enabled=%s",
                main_vcf,
                os.path.basename(main_vcf),
                has_csq_header,
                vep_enabled_from_pipeline,
            )

            # ── 2. FORMAT 문제 필드 정리 (동기 I/O → to_thread) ──
            cleaned_vcf = os.path.join(analysis_dir, f"{job.sample_name}.cleaned.vcf")
            main_mtime = os.path.getmtime(main_vcf)
            if (
                os.path.isfile(cleaned_vcf)
                and os.path.getsize(cleaned_vcf) > 0
                and os.path.getmtime(cleaned_vcf) >= main_mtime
            ):
                logger.info(
                    "  Reuse cleaned VCF (newer than main, skip FORMAT clean): %s",
                    cleaned_vcf,
                )
            else:
                logger.info(f"  Cleaning VCF FORMAT fields: {main_vcf} → {cleaned_vcf}")
                await asyncio.to_thread(clean_vcf_remove_formats, main_vcf, cleaned_vcf)
            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested after VCF clean")
                return False

            cleaned_mtime = os.path.getmtime(cleaned_vcf)

            # ── 3. snpEff annotation (선택적 – VEP VCF이면 건너뜀) ──
            annotated_vcf = cleaned_vcf
            snp_jar = settings.snpeff_jar
            if is_vep_vcf:
                # VEP annotation이 이미 완료된 경우 snpEff 실행 불필요
                logger.info("  VEP-annotated VCF detected: skipping snpEff step")
            elif snp_jar and os.path.exists(snp_jar):
                snpeff_vcf = os.path.join(analysis_dir, f"{job.sample_name}.snpeff.vcf")
                try:
                    if (
                        os.path.isfile(snpeff_vcf)
                        and os.path.getsize(snpeff_vcf) > 0
                        and os.path.getmtime(snpeff_vcf) >= cleaned_mtime
                    ):
                        annotated_vcf = snpeff_vcf
                        logger.info(
                            "  Reuse snpEff VCF (not older than cleaned VCF, skip snpEff): %s",
                            snpeff_vcf,
                        )
                    else:
                        logger.info("  Running snpEff annotation...")
                        await run_snpeff(
                            cleaned_vcf,
                            snpeff_vcf,
                            snp_jar,
                            settings.snpeff_db,
                            data_dir=settings.snpeff_data_dir,
                        )
                        annotated_vcf = snpeff_vcf
                        logger.info(f"  snpEff annotation complete: {snpeff_vcf}")
                except Exception as e:
                    logger.warning(f"  snpEff failed, continuing without: {e}")
            elif not is_vep_vcf and not snp_jar:
                logger.info("  snpEff not configured, skipping")
            elif not is_vep_vcf:
                logger.info("  snpEff JAR missing at %s, skipping", snp_jar)

            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested after snpEff step")
                return False

            # ── 4. BED 기반 필터링 준비 (job → CARRIER_DEFAULT_* → data/bed 탐색) ──
            backbone_bed_path = self._resolve_carrier_bed(
                job,
                "backbone_bed",
                "backbone",
                settings.carrier_default_backbone_bed,
            )
            backbone_regions = load_bed_regions(backbone_bed_path) if backbone_bed_path else {}
            logger.info(f"  Backbone BED: {backbone_bed_path or 'not configured'}")

            disease_bed_path = self._resolve_carrier_bed(
                job,
                "disease_bed",
                "disease",
                settings.carrier_default_disease_bed,
            )
            disease_regions = load_bed_regions(disease_bed_path) if disease_bed_path else {}
            logger.info(f"  Disease BED: {disease_bed_path or 'not configured'}")

            # ── 5. Annotator 초기화 (config 설정 활용) ──
            gnomad_dir = job.params.get("gnomad_dir") or (
                os.path.dirname(settings.gnomad_vcf) if settings.gnomad_vcf else
                (settings.gnomad_dir or "")
            )

            annotator = VariantAnnotator(
                clinvar_vcf=settings.clinvar_vcf or "",
                gnomad_dir=gnomad_dir,
                gnomad_genomes_glob=settings.gnomad_genomes_glob,
                gnomad_exomes_glob=settings.gnomad_exomes_glob,
                clingen_tsv=job.params.get("clingen_tsv") or settings.clingen_tsv or "",
                mane_gff=job.params.get("mane_gff") or settings.mane_gff or "",
                gene_bed=job.params.get("gene_bed") or settings.gene_bed or "",
                hpo_gene_file=settings.hpo_gene_file or "",
                curated_db=settings.curated_variants_db or "",
                hgmd_vcf=settings.hgmd_vcf or "",
                disease_gene_json=self._resolve_disease_gene_json_path(),
            )

            # ── 6. 필터 설정 구성 ──
            # HPO 유전자 필터 (job.params에서 HPO 텍스트가 제공된 경우)
            hpo_genes = set()
            hpo_text = job.params.get("hpo_terms", "")
            if hpo_text and annotator.hpo:
                from .annotator import HPOAnnotator
                hpo_ids = HPOAnnotator.normalize_hpo_terms(hpo_text)
                if hpo_ids:
                    hpo_genes = annotator.hpo.get_genes_for_hpo_list(hpo_ids)
                    logger.info(f"  HPO filter: {len(hpo_ids)} terms → {len(hpo_genes)} genes")

            # 유전자 필터 (job.params에서 유전자 목록이 제공된 경우)
            gene_filter_set = set()
            gene_filter_text = job.params.get("gene_filter", "")
            if gene_filter_text:
                gene_filter_set = {g.strip().upper() for g in gene_filter_text.split(",") if g.strip()}
                logger.info(f"  Gene filter: {len(gene_filter_set)} genes")

            # max AF 필터
            max_af = job.params.get("max_af")
            if max_af is not None:
                max_af = float(max_af)
                logger.info(f"  Max AF filter: {max_af}")

            # ClinVar 필터
            clinvar_filter_text = job.params.get("clinvar_filter", "")
            clinvar_filter = set()
            if clinvar_filter_text:
                clinvar_filter = {f.strip() for f in clinvar_filter_text.split(",") if f.strip()}
                logger.info(f"  ClinVar filter: {clinvar_filter}")

            restrict_disease = job.params.get("restrict_to_disease_bed")
            if restrict_disease is None:
                # Default False: panel scope = interpretation_genes post-filter (wes_panels), not per-panel BED.
                # Set True only if you intentionally want coordinate pre-filter from disease_bed_regions.
                restrict_disease = False
            filter_config = VariantFilterConfig(
                hpo_genes=hpo_genes,
                gene_filter_set=gene_filter_set,
                max_af=max_af,
                clinvar_filter=clinvar_filter,
                exclude_clinvar_conflicts=bool(job.params.get("exclude_clinvar_conflicts", False)),
                require_protein_altering=bool(job.params.get("require_protein_altering", True)),
                backbone_bed_regions=backbone_regions,
                disease_bed_regions=disease_regions,
                restrict_to_disease_bed=bool(restrict_disease),
            )

            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested before VCF parse")
                return False

            # ── 6b. VEP CSQ 사전 파싱 (VEP VCF인 경우) ──
            vep_annotations = None
            if is_vep_vcf:
                logger.info("  Pre-parsing VEP CSQ annotations from VCF...")
                try:
                    vep_annotations = await asyncio.to_thread(
                        extract_vep_annotations_from_vcf, annotated_vcf
                    )
                    logger.info(f"  VEP CSQ pre-parsing complete: {len(vep_annotations)} variants")
                except Exception as e:
                    logger.warning(f"  VEP CSQ pre-parsing failed, falling back to snpEff path: {e}")
                    vep_annotations = None

            # ── 7. 통합 VCF 파싱 파이프라인 실행 (CPU+I/O 블로킹 → to_thread) ──
            logger.info(
                "[process_results] parse_vcf_variants input (variants→result.json): %s | main_vcf was: %s",
                annotated_vcf,
                main_vcf,
            )
            logger.info(f"  Starting integrated VCF parsing pipeline...")
            annotated_variants, acmg_results, parse_stats = await asyncio.to_thread(
                parse_vcf_variants,
                vcf_path=annotated_vcf,
                annotator=annotator,
                filter_config=filter_config,
                acmg_classifier=None,
                vep_annotations=vep_annotations,
            )

            logger.info(
                f"  VCF parsing complete: "
                f"{parse_stats.get('total_records', 0)} records → "
                f"{parse_stats.get('final_count', 0)} variants after filtering "
                f"(disease/panel BED excluded: {parse_stats.get('disease_panel_filtered_count', 0)})"
            )

            if parse_stats.get("warnings"):
                for w in parse_stats["warnings"]:
                    logger.warning(f"  VCF parse warning: {w}")

            if parse_stats.get("total_records", 0) > 0 and parse_stats.get("final_count", 0) == 0:
                logger.warning(
                    "  VCF parse retained 0 variants (%s records, protein_altering prefilter=%s). "
                    "Often BED/HPO/gene_filter, max_af, ClinVar filter, or missing CSQ; check parse_stats.",
                    parse_stats.get("total_records"),
                    parse_stats.get("protein_altering_count"),
                )

            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested after VCF parse")
                return False

            # ── 8. ACMG 분류 ──────────────────────────────────────────
            # VUS 변이는 literature 검색 후 AI로 강화 (AI 활성화 시)
            from .acmg import classify_variant

            use_ai = settings.acmg_ai_enabled
            ai_provider = getattr(settings, "acmg_ai_provider", "gemini")
            if ai_provider == "gemini":
                acmg_api_key = getattr(settings, "gemini_api_key", "")
            else:
                acmg_api_key = getattr(settings, "acmg_ai_api_key", "")

            # literature 검색은 비동기로 사전 실행 (VUS 후보 기준 배치)
            literature_map: dict = {}
            if settings.literature_enabled and use_ai and annotated_variants:
                try:
                    from .literature import search_variant_literature
                    lit_tasks = []
                    lit_keys = []
                    for ann in annotated_variants:
                        gene = ann.get("gene", "")
                        hgvsc = ann.get("hgvsc", "")
                        hgvsp = ann.get("hgvsp", "")
                        if gene:
                            lit_tasks.append(
                                search_variant_literature(
                                    gene=gene,
                                    hgvsc=hgvsc,
                                    hgvsp=hgvsp,
                                    effect=ann.get("effect", ""),
                                    max_results=settings.literature_max_results,
                                    ncbi_email=settings.ncbi_email,
                                    ncbi_api_key=settings.ncbi_api_key,
                                    ncbi_tool=settings.ncbi_tool,
                                )
                            )
                            lit_keys.append(f"{gene}:{hgvsc or hgvsp}")
                    if lit_tasks:
                        logger.info(f"  Searching literature for {len(lit_tasks)} variants (semaphore=5)...")
                        # asyncio.gather 대신 세마포어로 동시성 제한 (NCBI 429 방지)
                        sem = asyncio.Semaphore(5)
                        async def _bounded(key, coro):
                            async with sem:
                                return key, await coro
                        bounded = [_bounded(k, t) for k, t in zip(lit_keys, lit_tasks)]
                        for fut in asyncio.as_completed(bounded):
                            try:
                                key, result = await fut
                                if isinstance(result, dict):
                                    literature_map[key] = result
                            except Exception as _lit_err:
                                logger.debug("[carrier_screening] Literature lookup failed: %s", _lit_err)
                        logger.info(
                            f"  Literature search complete: "
                            f"{sum(1 for v in literature_map.values() if v.get('total_found', 0) > 0)} "
                            f"variants with papers"
                        )
                except Exception as e:
                    logger.warning(f"  Literature search failed (non-fatal): {e}")

            if not acmg_results or all(r.get("final_classification") == "VUS" for r in acmg_results):
                logger.info(f"  Running ACMG classification for {len(annotated_variants)} variants...")
                acmg_results = []
                for ann in annotated_variants:
                    if _job_stop_requested(job.order_id):
                        logger.info("[process_results] User stop requested during ACMG step")
                        return False
                    gene = ann.get("gene", "")
                    hgvsc = ann.get("hgvsc", "")
                    hgvsp = ann.get("hgvsp", "")
                    lit = literature_map.get(f"{gene}:{hgvsc or hgvsp}")
                    acmg = await classify_variant(
                        ann,
                        use_ai=use_ai,
                        api_key=acmg_api_key,
                        provider=ai_provider,
                        literature=lit,
                    )
                    # 문헌 요약을 ACMG 결과에 포함 (result.json에서 참조 가능)
                    if lit:
                        acmg["literature"] = {
                            "total_found": lit.get("total_found", 0),
                            "from_cache": lit.get("from_cache", False),
                            "pmids": lit.get("pmids", []),
                            "top_titles": [
                                a.get("title", "")
                                for a in lit.get("articles", [])[:3]
                                if a.get("title")
                            ],
                        }
                    acmg_results.append(acmg)

            if should_apply_interpretation_post_filter(job):
                ig = interpretation_gene_set_for_job(job)
                if ig:
                    n_before = len(annotated_variants)
                    annotated_variants, acmg_results = _filter_variants_by_interpretation_genes(
                        annotated_variants, acmg_results, ig
                    )
                    logger.info(
                        "  Interpretation post-filter: %d gene(s) in set → %d variant(s) (was %d)",
                        len(ig),
                        len(annotated_variants),
                        n_before,
                    )

            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested before QC extract")
                return False

            # ── 9. QC 메트릭스 추출 (동기 I/O → to_thread) ──
            logger.info(f"  Extracting QC metrics...")
            qc_extras = self._qc_extra_search_dirs(job, analysis_dir)
            prior_ad = (job.params or {}).get("_prior_reuse_analysis_dir")
            if prior_ad:
                pad = str(prior_ad).strip()
                if pad and os.path.isdir(pad) and pad not in qc_extras:
                    qc_extras = [pad] + list(qc_extras)
            qc_sn = job.sample_name
            if job.params.get("_prior_reuse") and (job.params or {}).get(
                "_prior_reuse_sample_name"
            ):
                qc_sn = (
                    str((job.params or {}).get("_prior_reuse_sample_name") or "").strip()
                    or qc_sn
                )
            qc_more = self._qc_more_roots_from_outputs(job)
            qc_summary = await asyncio.to_thread(
                extract_qc_summary,
                analysis_dir,
                qc_sn,
                qc_extras,
                qc_more,
            )
            if job.params.get("_prior_reuse"):
                qc_summary["sample_name"] = job.sample_name

            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested after QC extract")
                return False

            # ── 10. Disease BED 정보 ──
            disease_bed_info = {}
            if disease_bed_path:
                total_genes = set()
                for regions in disease_regions.values():
                    for _, _, name in regions:
                        if name:
                            total_genes.add(name)
                disease_bed_info = {
                    "bed_file": os.path.basename(disease_bed_path),
                    "total_regions": sum(len(v) for v in disease_regions.values()),
                    "total_genes": len(total_genes),
                    "genes": sorted(total_genes),
                }

            # ── 11. 필터 요약 ──
            wes_pid = (job.params or {}).get("wes_panel_id")
            if not wes_pid and isinstance((job.params or {}).get("carrier"), dict):
                wes_pid = job.params["carrier"].get("wes_panel_id")
            wes_plabel = None
            if wes_pid:
                from ..wes_panels import get_panel_by_id

                _wp = get_panel_by_id(str(wes_pid))
                wes_plabel = (_wp.get("label") if _wp else None) or None

            _ig_set = interpretation_gene_set_for_job(job)
            _post_applied = should_apply_interpretation_post_filter(job) and bool(_ig_set)

            _prior_oid = (
                ((job.params or {}).get("carrier") or {}).get("prior_order_id")
                if job.params.get("_prior_reuse")
                else None
            )
            filter_summary = {
                "backbone_bed": os.path.basename(backbone_bed_path) if backbone_bed_path else None,
                "disease_bed": os.path.basename(disease_bed_path) if disease_bed_path else None,
                "restrict_to_disease_bed": bool(restrict_disease),
                "disease_panel_filtered_count": parse_stats.get(
                    "disease_panel_filtered_count", 0
                ),
                "hpo_terms_count": len(hpo_genes) if hpo_genes else 0,
                "gene_filter_count": len(gene_filter_set) if gene_filter_set else 0,
                "max_af": max_af,
                "clinvar_filter": list(clinvar_filter) if clinvar_filter else None,
                "require_protein_altering": filter_config.require_protein_altering,
                "wes_panel_id": str(wes_pid).strip() if wes_pid else None,
                "wes_panel_label": wes_plabel,
                "panel_filter_after_analysis": bool((job.params or {}).get("panel_filter_after_analysis")),
                "interpretation_post_filter_genes": len(_ig_set) if _post_applied else 0,
                "interpretation_post_filter_applied": _post_applied,
                "prior_order_id": str(_prior_oid).strip() if _prior_oid else None,
                "prior_pipeline_reuse": bool(job.params.get("_prior_reuse")),
                "include_apoe_pgx_requested": bool((job.params or {}).get("include_apoe_pgx")),
            }

            if _job_stop_requested(job.order_id):
                logger.info("[process_results] User stop requested before writing result.json")
                return False

            # ── 12. result.json 생성 (동기 I/O → to_thread) ──
            logger.info(f"  Generating result.json...")
            seen_mr = set()
            metric_image_roots: List[str] = []
            prior_out = (job.params or {}).get("_prior_reuse_output_dir")
            _mi_first: List[str] = []
            if prior_out:
                po = str(prior_out).strip()
                if po and os.path.isdir(po):
                    _mi_first.append(po)
            for p in _mi_first + [artifact_output_dir, portal_output_dir] + list(qc_more or []):
                if not p or not os.path.isdir(p):
                    continue
                try:
                    rp = os.path.realpath(p)
                except OSError:
                    continue
                if rp in seen_mr:
                    continue
                seen_mr.add(rp)
                metric_image_roots.append(p)
            dark_roots = self._qc_extra_search_dirs(job, analysis_dir)
            _pgx_allow = pgx_portal_gene_allowlist_for_job(job)
            result_json_path = await asyncio.to_thread(
                generate_result_json,
                annotated_variants=annotated_variants,
                acmg_results=acmg_results,
                qc_summary=qc_summary,
                sample_name=job.sample_name,
                order_id=job.order_id,
                disease_bed_info=disease_bed_info,
                filter_summary=filter_summary,
                parse_stats=parse_stats,
                output_dir=portal_output_dir,
                metric_image_search_roots=metric_image_roots,
                analysis_dir=analysis_dir,
                dark_genes_extra_roots=dark_roots,
                review_build_metadata={
                    "main_vcf": main_vcf,
                    "annotated_vcf": annotated_vcf,
                    "vep_preparsed_keys": len(vep_annotations) if vep_annotations else 0,
                    "vep_annotation_enabled": bool(is_vep_vcf),
                },
                pgx_panel_gene_symbols=_pgx_allow,
            )

            # ── 13. variants.tsv 생성 ──
            logger.info(f"  Generating variants.tsv...")

            def _write_tsv():
                with open(result_json_path, "r") as f:
                    result_data = json.load(f)
                generate_variants_tsv(result_data.get("variants", []), portal_output_dir)

            await asyncio.to_thread(_write_tsv)

            # ── 14. QC summary JSON 저장 ──
            def _write_qc():
                qc_path = os.path.join(portal_output_dir, "qc_summary.json")
                with open(qc_path, "w", encoding="utf-8") as f:
                    json.dump(qc_summary, f, ensure_ascii=False, indent=2, default=str)

            await asyncio.to_thread(_write_qc)

            logger.info(
                f"[process_results] Complete: portal result → {portal_output_dir} "
                f"({len(annotated_variants)} variants, "
                f"{len(acmg_results)} ACMG classifications)"
            )
            return True

        except ImportError as e:
            logger.error(f"Missing dependency for annotation: {e}")
            return False
        except Exception as e:
            logger.error(f"Result processing failed: {e}", exc_info=True)
            return False

    # ─── Step 5: 출력 파일 목록 ────────────────────────────

    async def get_output_files(self, job: Job) -> List[OutputFile]:
        """Portal에 업로드할 파일 목록을 반환합니다."""
        dirs = self._get_dirs(job)
        review_out = carrier_report_output_dir(job)
        files = []

        # result.json (필수)
        result_json = os.path.join(review_out, "result.json")
        if os.path.exists(result_json):
            files.append(OutputFile(
                file_path=result_json,
                file_type="review_json",
                file_name="result.json",
                content_type="application/json",
            ))

        # qc_summary.json
        qc_json = os.path.join(review_out, "qc_summary.json")
        if os.path.exists(qc_json):
            files.append(OutputFile(
                file_path=qc_json,
                file_type="qc_json",
                file_name="qc_summary.json",
                content_type="application/json",
            ))

        # variants.tsv
        variants_tsv = os.path.join(review_out, "variants.tsv")
        if os.path.exists(variants_tsv):
            files.append(OutputFile(
                file_path=variants_tsv,
                file_type="variants_tsv",
                file_name="variants.tsv",
                content_type="text/tab-separated-values",
            ))

        # IGV 스냅샷 이미지
        snapshot_dirs = [
            os.path.join(review_out, "snapshots"),
            os.path.join(dirs["analysis"], "snapshots"),
        ]
        for snap_dir in snapshot_dirs:
            if os.path.isdir(snap_dir):
                for img in glob.glob(os.path.join(snap_dir, "*.png")):
                    files.append(OutputFile(
                        file_path=img,
                        file_type="igv_snapshot",
                        file_name=os.path.basename(img),
                        content_type="image/png",
                    ))
                break  # 첫 번째 존재하는 디렉토리만 사용

        logger.info(f"Output files for upload: {len(files)} files")
        return files

    # ─── Step 6: 리포트 생성 (리뷰어 확정 후) ──────────────

    def _generate_report_sync(
        self,
        job: Job,
        confirmed_variants: List[Dict],
        reviewer_info: Dict,
        patient_info: Optional[Dict] = None,
        partner_info: Optional[Dict] = None,
        languages: Optional[List[str]] = None,
    ) -> bool:
        """WeasyPrint 등 동기 작업 — asyncio.to_thread 에서 실행."""
        try:
            from .report import (
                CARRIER_PDF_SOLO_KINDS,
                generate_report_json,
                generate_report_pdf,
                carrier_report_template_kind,
                report_languages_from_order,
                _carrier_order_flat,
            )
            from .review import extract_qc_summary

            output_dir, analysis_dir = carrier_report_generation_paths(job)

            # QC 요약
            qc_extras = self._qc_extra_search_dirs(job, analysis_dir)
            qc_more = self._qc_more_roots_from_outputs(job)
            qc_summary = extract_qc_summary(
                analysis_dir, job.sample_name, qc_extras, qc_more
            )

            raw_params = job.params or {}
            kind = carrier_report_template_kind(raw_params)
            langs = report_languages_from_order(raw_params)
            params = _carrier_order_flat(raw_params)
            if kind is None or langs is None:
                logger.warning(
                    f"[generate_report] Unsupported carrier PDF order: {job.order_id} "
                    f"kind={kind} langs={langs} params_keys={list(params.keys())}"
                )
                return False

            languages = langs

            # Solo PDFs: carrier / exome / proactive / pgx — 파트너 섹션 없음
            # CouplesCarrier: carrier_couples_<lang>.html — patient2_* 또는 폼 partner
            if kind in CARRIER_PDF_SOLO_KINDS:
                partner_info = None
            else:
                merged: Dict[str, Any] = {}
                if partner_info and (partner_info.get("name") or "").strip():
                    merged = dict(partner_info)
                else:
                    p2n = (params.get("patient2_name") or "").strip()
                    if p2n:
                        merged["name"] = p2n
                    if params.get("patient2_birth"):
                        merged["dob"] = params["patient2_birth"]
                    if params.get("patient2_gender"):
                        merged["gender"] = params["patient2_gender"]
                if not (merged.get("name") or "").strip():
                    logger.warning(
                        f"[generate_report] Couples report missing partner name: {job.order_id}"
                    )
                    return False
                partner_info = merged

            pi = dict(patient_info) if patient_info else {}
            if not (pi.get("name") or "").strip():
                pi["name"] = (params.get("patient_name") or "").strip() or job.sample_name
            if not pi.get("dob") and params.get("patient_birth"):
                pi["dob"] = params["patient_birth"]
            if not pi.get("gender") and params.get("patient_gender"):
                pi["gender"] = params["patient_gender"]
            patient_info = pi

            # 템플릿: 명시 → output_root/report_templates → data/report_templates → REPORT_TEMPLATE_DIR → 파이프라인
            template_dir = resolve_carrier_pdf_template_dir()

            logger.info(
                f"[generate_report] Generating report for {job.order_id} "
                f"(languages={languages}, template_dir={template_dir})"
            )

            # report.json 생성
            report_json_path = generate_report_json(
                order_id=job.order_id,
                sample_name=job.sample_name,
                confirmed_variants=confirmed_variants,
                reviewer_info=reviewer_info,
                qc_summary=qc_summary,
                output_dir=output_dir,
                patient_info=patient_info,
                partner_info=partner_info,
                order_params=job.params,
                report_language=(languages[0] if languages else "EN"),
                disease_gene_json=self._resolve_disease_gene_json_path(),
                gene_knowledge_db=settings.gene_knowledge_db or None,
                gene_knowledge_enrich_on_report=bool(
                    getattr(settings, "gene_knowledge_enrich_on_report", True)
                ),
                gene_knowledge_gemini_on_report=bool(
                    getattr(settings, "gene_knowledge_gemini_on_report", True)
                ),
                gemini_api_key=settings.gemini_api_key or None,
                gene_knowledge_gemini_model=getattr(
                    settings, "gene_knowledge_gemini_model", "gemini-2.5-flash"
                ),
                extra_result_json_paths=_extra_result_json_paths_for_carrier_report(job),
                pdf_template_kind=kind,
            )
            logger.info(f"  Generated report.json: {report_json_path}")

            # 다국어 PDF 리포트 생성
            pdf_paths = generate_report_pdf(
                report_json_path=report_json_path,
                output_dir=output_dir,
                template_dir=template_dir,
                languages=languages,
                extra_result_json_paths=_extra_result_json_paths_for_carrier_report(job),
            )

            for pdf_path in pdf_paths:
                logger.info(f"  Generated PDF: {pdf_path}")

            logger.info(
                f"[generate_report] Complete: {output_dir} "
                f"({len(pdf_paths)} PDF files generated)"
            )
            return True

        except Exception as e:
            logger.error(f"Report generation failed: {e}", exc_info=True)
            return False

    async def generate_report(
        self,
        job: Job,
        confirmed_variants: List[Dict],
        reviewer_info: Dict,
        patient_info: Optional[Dict] = None,
        partner_info: Optional[Dict] = None,
        languages: Optional[List[str]] = None,
    ) -> bool:
        """
        리뷰어 확정 후 최종 리포트를 생성합니다.
        service-daemon의 /order/{order_id}/report 엔드포인트에서 호출됩니다.

        WeasyPrint PDF 생성은 이벤트 루프를 막지 않도록 스레드에서 실행합니다.
        """
        return await asyncio.to_thread(
            self._generate_report_sync,
            job,
            confirmed_variants,
            reviewer_info,
            patient_info,
            partner_info,
            languages,
        )

    # ─── Lifecycle Hooks ───────────────────────────────────

    async def on_job_start(self, job: Job):
        """작업 시작 시 호출"""
        logger.info(f"[carrier_screening] Job started: {job.order_id} (sample: {job.sample_name})")

    async def on_job_complete(self, job: Job):
        """작업 완료 시 호출"""
        logger.info(f"[carrier_screening] Job completed: {job.order_id}")

    async def on_job_failed(self, job: Job, error: str):
        """작업 실패 시 호출"""
        logger.error(f"[carrier_screening] Job failed: {job.order_id} - {error}")

    def sync_is_complete(self, job: Job) -> bool:
        """
        daemon 재시작 복구 시 파이프라인이 실제로 완료됐는지 동기적으로 확인.
        result.json 이 존재하면 process_results 까지 완료된 것으로 판정.
        """
        rj = carrier_result_json_path(job)
        if rj and os.path.isfile(rj):
            logger.info(
                "[carrier_screening] sync_is_complete: found result.json at %s → marking COMPLETED", rj
            )
            return True
        return False

    async def cleanup(self, job: Job):
        """작업 정리"""
        # Nextflow work 디렉토리 정리 (선택적)
        dirs = self._get_dirs(job)
        work_dir = os.path.join(dirs["analysis"], "work")
        if os.path.isdir(work_dir) and job.params.get("cleanup_work", False):
            import shutil
            try:
                shutil.rmtree(work_dir)
                logger.info(f"Cleaned up Nextflow work dir: {work_dir}")
            except Exception as e:
                logger.warning(f"Failed to cleanup work dir: {e}")
