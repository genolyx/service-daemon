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
from typing import Dict, Any, List, Optional, Tuple

from ..base import ServicePlugin
from ...config import settings
from ...models import Job, OutputFile

logger = logging.getLogger(__name__)


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

    def validate_params(self, params: Dict[str, Any]) -> Tuple[bool, str]:
        """서비스별 파라미터 유효성 검사"""
        # 필수 파라미터 없음 (기본값 사용 가능)
        return True, ""

    # ─── 디렉토리 구조 ─────────────────────────────────────

    def _get_dirs(self, job: Job) -> Dict[str, str]:
        """
        carrier-screening 디렉토리 구조:
            carrier-screening/fastq/<work_dir>/<sample_name>/
            carrier-screening/analysis/<work_dir>/<sample_name>/
            carrier-screening/output/<work_dir>/<sample_name>/
            carrier-screening/log/<work_dir>/<sample_name>/
        """
        base = os.path.join(settings.base_dir, "carrier-screening")
        return {
            "fastq": job.fastq_dir or os.path.join(base, "fastq", job.work_dir, job.sample_name),
            "analysis": job.analysis_dir or os.path.join(base, "analysis", job.work_dir, job.sample_name),
            "output": job.output_dir or os.path.join(base, "output", job.work_dir, job.sample_name),
            "log": job.log_dir or os.path.join(base, "log", job.work_dir, job.sample_name),
        }

    # ─── Step 1: 입력 준비 ─────────────────────────────────

    async def prepare_inputs(self, job: Job) -> bool:
        """FASTQ 파일 확인 및 디렉토리 구조 생성"""
        dirs = self._get_dirs(job)

        # 디렉토리 생성
        for d in dirs.values():
            os.makedirs(d, exist_ok=True)

        # Job 경로 업데이트
        job.fastq_dir = dirs["fastq"]
        job.analysis_dir = dirs["analysis"]
        job.output_dir = dirs["output"]
        job.log_dir = dirs["log"]

        # FASTQ 파일 확인
        r1_path, r2_path = self._find_fastq_files(dirs["fastq"], job.sample_name)

        if r1_path and r2_path:
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
                job.fastq_r1_path = r1_dest
                job.fastq_r2_path = r2_dest
                return True
            else:
                logger.error("Failed to download FASTQ files")
                return False

        # 지정된 로컬 경로 확인
        if job.fastq_r1_path and job.fastq_r2_path:
            if os.path.exists(job.fastq_r1_path) and os.path.exists(job.fastq_r2_path):
                return True

        logger.error(f"No FASTQ files found for sample {job.sample_name}")
        return False

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
            cmd = f"wget -q -O '{dest}' '{url}'"
            proc = await asyncio.create_subprocess_shell(
                cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
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

    async def get_pipeline_command(self, job: Job) -> str:
        """Nextflow 파이프라인 실행 명령을 생성합니다."""
        pipeline_dir = settings.carrier_screening_pipeline_dir
        main_nf = os.path.join(pipeline_dir, "bin", "main.nf")
        nf_config = os.path.join(pipeline_dir, "bin", "nextflow.config")

        # Nextflow 실행 명령 구성
        cmd_parts = [
            settings.nextflow_executable,
            "run", main_nf,
        ]

        # config 파일
        if os.path.exists(nf_config):
            cmd_parts.extend(["-c", nf_config])

        # 파라미터
        params = {
            "sample_name": job.sample_name,
            "fastq_r1": job.fastq_r1_path,
            "fastq_r2": job.fastq_r2_path,
            "outdir": job.analysis_dir,
        }

        # Reference 파일 (설정에서 가져오기)
        if settings.ref_fasta:
            params["ref"] = settings.ref_fasta
        if settings.ref_bwa_indices:
            params["bwa_index"] = settings.ref_bwa_indices

        # BED 파일 (서비스 파라미터에서 가져오기)
        bed_dir = os.path.join(pipeline_dir, "data", "bed")
        backbone_bed = job.params.get("backbone_bed") or self._find_default_bed(bed_dir, "backbone")
        if backbone_bed:
            params["backbone_bed"] = backbone_bed

        # 추가 파라미터 (job.params에서 오버라이드)
        for key in ("pon_tar", "target_bed", "disease_bed", "cnv_bed"):
            if job.params.get(key):
                params[key] = job.params[key]

        # 파라미터를 명령줄에 추가
        for key, val in params.items():
            if val:
                cmd_parts.append(f"--{key}")
                cmd_parts.append(str(val))

        # 작업 디렉토리
        work_dir = os.path.join(job.analysis_dir, "work")
        cmd_parts.extend(["-work-dir", work_dir])

        # 리포트
        report_path = os.path.join(job.log_dir, "nextflow_report.html")
        cmd_parts.extend(["-with-report", report_path])

        return " ".join(cmd_parts)

    def _find_default_bed(self, bed_dir: str, prefix: str) -> Optional[str]:
        """기본 BED 파일을 찾습니다."""
        if not os.path.isdir(bed_dir):
            return None
        candidates = glob.glob(os.path.join(bed_dir, f"{prefix}*.bed"))
        if candidates:
            return candidates[0]
        return None

    # ─── Step 3: 완료 확인 ─────────────────────────────────

    async def check_completion(self, job: Job) -> bool:
        """파이프라인 완료를 확인합니다 (VCF 파일 존재 여부)."""
        dirs = self._get_dirs(job)
        analysis_dir = dirs["analysis"]

        # VCF 파일 찾기
        vcf_files = glob.glob(os.path.join(analysis_dir, "**", "*.vcf.gz"), recursive=True)
        vcf_files += glob.glob(os.path.join(analysis_dir, "**", "*.vcf"), recursive=True)

        if not vcf_files:
            logger.error(f"No VCF files found in {analysis_dir}")
            return False

        # 주요 VCF 파일 확인 (genotyped, filtered)
        main_vcf = None
        for pattern in [
            f"*{job.sample_name}*filtered*.vcf*",
            f"*{job.sample_name}*genotyped*.vcf*",
            f"*{job.sample_name}*.vcf*",
        ]:
            matches = glob.glob(os.path.join(analysis_dir, "**", pattern), recursive=True)
            if matches:
                main_vcf = matches[0]
                break

        if main_vcf:
            logger.info(f"Pipeline completed, main VCF: {main_vcf}")
            job.params["main_vcf"] = main_vcf
            return True

        logger.warning(f"VCF files found but no main VCF identified: {vcf_files}")
        job.params["main_vcf"] = vcf_files[0]
        return True

    # ─── Step 4: 결과 후처리 (Annotation) ──────────────────

    async def process_results(self, job: Job) -> bool:
        """
        VCF → BED 필터 → Annotation → ACMG 분류 → result.json 생성

        통합 parse_vcf_variants 파이프라인을 사용하여 처리합니다.

        전체 annotation 파이프라인:
            1. VCF FORMAT 문제 필드 정리
            2. snpEff annotation (선택적)
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

            dirs = self._get_dirs(job)
            analysis_dir = dirs["analysis"]
            output_dir = dirs["output"]
            os.makedirs(output_dir, exist_ok=True)

            # ── 1. VCF 파일 확인 ──
            main_vcf = job.params.get("main_vcf")
            if not main_vcf or not os.path.exists(main_vcf):
                logger.error(f"Main VCF not found: {main_vcf}")
                return False

            logger.info(f"[process_results] Starting annotation for {job.sample_name}")
            logger.info(f"  Main VCF: {main_vcf}")

            # ── 2. FORMAT 문제 필드 정리 ──
            cleaned_vcf = os.path.join(analysis_dir, f"{job.sample_name}.cleaned.vcf")
            logger.info(f"  Cleaning VCF FORMAT fields: {main_vcf} → {cleaned_vcf}")
            clean_vcf_remove_formats(main_vcf, cleaned_vcf)

            # ── 3. snpEff annotation (선택적) ──
            annotated_vcf = cleaned_vcf
            if settings.snpeff_jar and os.path.exists(settings.snpeff_jar):
                snpeff_vcf = os.path.join(analysis_dir, f"{job.sample_name}.snpeff.vcf")
                try:
                    logger.info(f"  Running snpEff annotation...")
                    run_snpeff(cleaned_vcf, snpeff_vcf, settings.snpeff_jar, settings.snpeff_db)
                    annotated_vcf = snpeff_vcf
                    logger.info(f"  snpEff annotation complete: {snpeff_vcf}")
                except Exception as e:
                    logger.warning(f"  snpEff failed, continuing without: {e}")
            else:
                logger.info("  snpEff not configured, skipping")

            # ── 4. BED 기반 필터링 준비 ──
            pipeline_dir = settings.carrier_screening_pipeline_dir
            bed_dir = os.path.join(pipeline_dir, "data", "bed")

            # backbone BED (전체 타겟 영역)
            backbone_bed_path = job.params.get("backbone_bed") or self._find_default_bed(bed_dir, "backbone")
            backbone_regions = load_bed_regions(backbone_bed_path) if backbone_bed_path else {}
            logger.info(f"  Backbone BED: {backbone_bed_path or 'not configured'}")

            # disease BED (질환 관련 유전자)
            disease_bed_path = job.params.get("disease_bed") or self._find_default_bed(bed_dir, "disease")
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
                disease_gene_json=settings.disease_gene_json or "",
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

            filter_config = VariantFilterConfig(
                hpo_genes=hpo_genes,
                gene_filter_set=gene_filter_set,
                max_af=max_af,
                clinvar_filter=clinvar_filter,
                exclude_clinvar_conflicts=bool(job.params.get("exclude_clinvar_conflicts", False)),
                require_protein_altering=bool(job.params.get("require_protein_altering", True)),
                backbone_bed_regions=backbone_regions,
                disease_bed_regions=disease_regions,
            )

            # ── 7. 통합 VCF 파싱 파이프라인 실행 ──
            logger.info(f"  Starting integrated VCF parsing pipeline...")
            annotated_variants, acmg_results, parse_stats = parse_vcf_variants(
                vcf_path=annotated_vcf,
                annotator=annotator,
                filter_config=filter_config,
                acmg_classifier=None,  # rule-based ACMG는 parse_vcf_variants 내부에서 처리
            )

            logger.info(
                f"  VCF parsing complete: "
                f"{parse_stats.get('total_records', 0)} records → "
                f"{parse_stats.get('final_count', 0)} variants after filtering"
            )

            if parse_stats.get("warnings"):
                for w in parse_stats["warnings"]:
                    logger.warning(f"  VCF parse warning: {w}")

            # ── 8. ACMG 분류 (parse_vcf_variants에서 처리되지 않은 경우) ──
            if not acmg_results or all(r.get("final_classification") == "VUS" for r in acmg_results):
                from .acmg import classify_variant
                logger.info(f"  Running ACMG classification for {len(annotated_variants)} variants...")
                acmg_results = []
                for ann in annotated_variants:
                    acmg = await classify_variant(ann, use_ai=False)
                    acmg_results.append(acmg)

            # ── 9. QC 메트릭스 추출 ──
            logger.info(f"  Extracting QC metrics...")
            qc_summary = extract_qc_summary(analysis_dir, job.sample_name)

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
            filter_summary = {
                "backbone_bed": os.path.basename(backbone_bed_path) if backbone_bed_path else None,
                "disease_bed": os.path.basename(disease_bed_path) if disease_bed_path else None,
                "hpo_terms_count": len(hpo_genes) if hpo_genes else 0,
                "gene_filter_count": len(gene_filter_set) if gene_filter_set else 0,
                "max_af": max_af,
                "clinvar_filter": list(clinvar_filter) if clinvar_filter else None,
                "require_protein_altering": filter_config.require_protein_altering,
            }

            # ── 12. result.json 생성 ──
            logger.info(f"  Generating result.json...")
            result_json_path = generate_result_json(
                annotated_variants=annotated_variants,
                acmg_results=acmg_results,
                qc_summary=qc_summary,
                sample_name=job.sample_name,
                order_id=job.order_id,
                disease_bed_info=disease_bed_info,
                filter_summary=filter_summary,
                parse_stats=parse_stats,
                output_dir=output_dir,
            )

            # ── 13. variants.tsv 생성 ──
            logger.info(f"  Generating variants.tsv...")
            with open(result_json_path, "r") as f:
                result_data = json.load(f)
            generate_variants_tsv(result_data.get("variants", []), output_dir)

            # ── 14. QC summary JSON 저장 ──
            qc_path = os.path.join(output_dir, "qc_summary.json")
            with open(qc_path, "w", encoding="utf-8") as f:
                json.dump(qc_summary, f, ensure_ascii=False, indent=2, default=str)

            logger.info(
                f"[process_results] Complete: {output_dir} "
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
        output_dir = dirs["output"]
        files = []

        # result.json (필수)
        result_json = os.path.join(output_dir, "result.json")
        if os.path.exists(result_json):
            files.append(OutputFile(
                file_path=result_json,
                file_type="review_json",
                file_name="result.json",
                content_type="application/json",
            ))

        # qc_summary.json
        qc_json = os.path.join(output_dir, "qc_summary.json")
        if os.path.exists(qc_json):
            files.append(OutputFile(
                file_path=qc_json,
                file_type="qc_json",
                file_name="qc_summary.json",
                content_type="application/json",
            ))

        # variants.tsv
        variants_tsv = os.path.join(output_dir, "variants.tsv")
        if os.path.exists(variants_tsv):
            files.append(OutputFile(
                file_path=variants_tsv,
                file_type="variants_tsv",
                file_name="variants.tsv",
                content_type="text/tab-separated-values",
            ))

        # IGV 스냅샷 이미지
        snapshot_dirs = [
            os.path.join(output_dir, "snapshots"),
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

        Args:
            job: 작업 정보
            confirmed_variants: 리뷰어가 확정한 변이 목록
            reviewer_info: 리뷰어 정보
            patient_info: 환자 정보 (선택)
            partner_info: 파트너 정보 (couple 검사 시, 선택)
            languages: 리포트 생성 언어 목록 (기본: config의 report_language_list)
        """
        try:
            from .report import generate_report_json, generate_report_pdf
            from .review import extract_qc_summary

            dirs = self._get_dirs(job)
            output_dir = dirs["output"]
            analysis_dir = dirs["analysis"]

            # QC 요약
            qc_summary = extract_qc_summary(analysis_dir, job.sample_name)

            # 리포트 언어 결정
            if languages is None:
                languages = settings.report_language_list

            # 템플릿 디렉토리 결정
            template_dir = settings.report_template_dir
            if not template_dir:
                # 파이프라인 내 templates 디렉토리 확인
                pipeline_template_dir = os.path.join(
                    settings.carrier_screening_pipeline_dir, "data", "templates"
                )
                if os.path.isdir(pipeline_template_dir):
                    template_dir = pipeline_template_dir

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
            )
            logger.info(f"  Generated report.json: {report_json_path}")

            # 다국어 PDF 리포트 생성
            pdf_paths = generate_report_pdf(
                report_json_path=report_json_path,
                output_dir=output_dir,
                template_dir=template_dir,
                languages=languages,
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
