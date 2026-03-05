"""
Carrier Screening Service Plugin

ServicePlugin 인터페이스 구현체.
전체 워크플로우를 조율합니다:

    1. prepare_inputs     : FASTQ 파일 확인/다운로드, 디렉토리 구조 생성
    2. get_pipeline_command: Nextflow 파이프라인 실행 명령 생성
    3. check_completion   : 파이프라인 완료 확인 (VCF 존재 여부)
    4. process_results    : VCF → BED 필터 → Annotation → ACMG → result.json 생성
    5. get_output_files   : Portal에 업로드할 파일 목록 반환
    6. generate_report    : 리뷰어 확정 후 report.json + PDF 생성
"""

import os
import re
import glob
import asyncio
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

        전체 annotation 파이프라인:
            1. VCF FORMAT 문제 필드 정리
            2. snpEff annotation (선택적)
            3. BED 기반 변이 필터링
            4. pysam으로 변이 파싱
            5. ClinVar, gnomAD, ClinGen, MANE annotation
            6. ACMG 분류
            7. QC 메트릭스 추출
            8. result.json 생성
        """
        try:
            from .vcf_parser import (
                clean_vcf_remove_formats, run_snpeff,
                get_annotation_layout, extract_variant_info,
                get_sample_metrics, is_protein_altering, is_nonref,
                load_bed_regions, variant_in_bed,
            )
            from .annotator import VariantAnnotator
            from .acmg import classify_variant
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

            # ── 2. FORMAT 문제 필드 정리 ──
            cleaned_vcf = os.path.join(analysis_dir, f"{job.sample_name}.cleaned.vcf")
            logger.info(f"Cleaning VCF FORMAT fields: {main_vcf} → {cleaned_vcf}")
            clean_vcf_remove_formats(main_vcf, cleaned_vcf)

            # ── 3. snpEff annotation (선택적) ──
            annotated_vcf = cleaned_vcf
            if settings.snpeff_jar and os.path.exists(settings.snpeff_jar):
                snpeff_vcf = os.path.join(analysis_dir, f"{job.sample_name}.snpeff.vcf")
                try:
                    run_snpeff(cleaned_vcf, snpeff_vcf, settings.snpeff_jar, settings.snpeff_db)
                    annotated_vcf = snpeff_vcf
                except Exception as e:
                    logger.warning(f"snpEff failed, continuing without: {e}")

            # ── 4. BED 기반 필터링 ──
            pipeline_dir = settings.carrier_screening_pipeline_dir
            bed_dir = os.path.join(pipeline_dir, "data", "bed")

            # backbone BED (전체 타겟 영역)
            backbone_bed_path = job.params.get("backbone_bed") or self._find_default_bed(bed_dir, "backbone")
            backbone_regions = load_bed_regions(backbone_bed_path) if backbone_bed_path else {}

            # disease BED (질환 관련 유전자)
            disease_bed_path = job.params.get("disease_bed") or self._find_default_bed(bed_dir, "disease")
            disease_regions = load_bed_regions(disease_bed_path) if disease_bed_path else {}

            # ── 5. Annotator 초기화 ──
            gnomad_dir = job.params.get("gnomad_dir") or os.path.dirname(settings.gnomad_vcf or "")
            clingen_tsv = job.params.get("clingen_tsv", "")
            mane_gff = job.params.get("mane_gff", "")
            gene_bed = job.params.get("gene_bed", "")

            annotator = VariantAnnotator(
                clinvar_vcf=settings.clinvar_vcf or "",
                gnomad_dir=gnomad_dir,
                clingen_tsv=clingen_tsv,
                mane_gff=mane_gff,
                gene_bed=gene_bed,
            )

            # ── 6. VCF 파싱 및 Annotation ──
            import pysam
            vcf = pysam.VariantFile(annotated_vcf)
            tag, gi, ti, ci, pi, ei = get_annotation_layout(vcf)

            annotated_variants = []
            acmg_results = []

            for rec in vcf:
                # non-reference 필터
                if not is_nonref(rec):
                    continue

                # BED 영역 필터 (backbone)
                if backbone_regions and not variant_in_bed(rec.chrom, rec.pos, backbone_regions):
                    continue

                for alt in (rec.alts or []):
                    alt = str(alt)

                    # 유전자/전사체/HGVS 추출
                    gene, transcript, hgvsc, hgvsp, effect = extract_variant_info(
                        rec, tag, gi, ti, ci, pi, ei,
                        gene_interval_lookup=annotator.gene_intervals.lookup if annotator.gene_intervals else None,
                    )

                    # 샘플 메트릭스
                    metrics = get_sample_metrics(rec, alt)

                    # 통합 annotation
                    ann = annotator.annotate(
                        chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=alt,
                        gene=gene, transcript=transcript,
                        hgvsc=hgvsc, hgvsp=hgvsp, effect=effect,
                        sample_metrics=metrics,
                    )

                    # ACMG 분류
                    acmg = await classify_variant(ann, use_ai=False)

                    annotated_variants.append(ann)
                    acmg_results.append(acmg)

            vcf.close()
            logger.info(f"Annotated {len(annotated_variants)} variants")

            # ── 7. QC 메트릭스 추출 ──
            qc_summary = extract_qc_summary(analysis_dir, job.sample_name)

            # ── 8. Disease BED 정보 ──
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

            # ── 9. result.json 생성 ──
            result_json_path = generate_result_json(
                annotated_variants=annotated_variants,
                acmg_results=acmg_results,
                qc_summary=qc_summary,
                sample_name=job.sample_name,
                order_id=job.order_id,
                disease_bed_info=disease_bed_info,
                output_dir=output_dir,
            )

            # ── 10. variants.tsv 생성 ──
            # result.json의 variants를 TSV로도 내보내기
            import json
            with open(result_json_path, "r") as f:
                result_data = json.load(f)
            generate_variants_tsv(result_data.get("variants", []), output_dir)

            # ── 11. QC summary JSON 저장 ──
            qc_path = os.path.join(output_dir, "qc_summary.json")
            with open(qc_path, "w", encoding="utf-8") as f:
                json.dump(qc_summary, f, ensure_ascii=False, indent=2, default=str)

            logger.info(f"Result processing complete: {output_dir}")
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

    async def generate_report(self, job: Job, confirmed_variants: List[Dict], reviewer_info: Dict) -> bool:
        """
        리뷰어 확정 후 최종 리포트를 생성합니다.
        service-daemon의 /order/{order_id}/report 엔드포인트에서 호출됩니다.
        """
        try:
            from .report import generate_report_json, generate_report_pdf
            from .review import extract_qc_summary

            dirs = self._get_dirs(job)
            output_dir = dirs["output"]
            analysis_dir = dirs["analysis"]

            # QC 요약
            qc_summary = extract_qc_summary(analysis_dir, job.sample_name)

            # report.json 생성
            report_json_path = generate_report_json(
                order_id=job.order_id,
                sample_name=job.sample_name,
                confirmed_variants=confirmed_variants,
                reviewer_info=reviewer_info,
                qc_summary=qc_summary,
                output_dir=output_dir,
            )

            # PDF 리포트 생성
            pdf_path = generate_report_pdf(report_json_path, output_dir)

            logger.info(f"Report generation complete: {output_dir}")
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
