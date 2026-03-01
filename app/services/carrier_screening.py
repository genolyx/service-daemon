"""
Carrier Screening Service Plugin

Carrier Screening 파이프라인 실행, Annotation, 리뷰 페이지 생성을 담당합니다.

워크플로우:
1. FASTQ 입력 준비
2. Nextflow 파이프라인 실행 (carrier-screening/bin/main.nf)
3. 완료 확인 (summary_report 존재 여부)
4. 결과 후처리:
   a. VCF Annotation (ClinVar, gnomAD, dbSNP)
   b. 리뷰 페이지용 JSON 생성
   c. 스냅샷 이미지 정리
5. 출력 파일 목록 반환 (Platform 업로드용)
"""

import os
import re
import json
import glob
import gzip
import shutil
import logging
import subprocess
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path

from .base import ServicePlugin
from app.config import settings
from app.models import Job, OutputFile

logger = logging.getLogger(__name__)


# ─── Annotation Helpers ────────────────────────────────────

class VCFAnnotator:
    """
    VCF 파일에 ClinVar, gnomAD, dbSNP 등의 annotation을 추가합니다.
    carrier-screening/Genomics_Pipeline/phenotype_portal/main.py의 로직을 재사용합니다.
    """

    def __init__(self):
        self.clinvar_vcf = settings.clinvar_vcf
        self.gnomad_vcf = settings.gnomad_vcf
        self.dbsnp_vcf = settings.dbsnp_vcf
        self.snpeff_jar = settings.snpeff_jar
        self.snpeff_db = settings.snpeff_db

    def annotate_with_snpeff(self, input_vcf: str, output_vcf: str) -> bool:
        """snpEff를 사용한 기능적 annotation"""
        if not self.snpeff_jar or not os.path.exists(self.snpeff_jar):
            logger.warning("snpEff JAR not found, skipping functional annotation")
            return True  # 선택적이므로 실패하지 않음

        try:
            cmd = [
                "java", "-Xmx4g", "-jar", self.snpeff_jar,
                "ann", "-v", self.snpeff_db,
                input_vcf
            ]
            with open(output_vcf, "w") as out_f:
                subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE,
                               check=True, timeout=1800)
            logger.info(f"snpEff annotation complete: {output_vcf}")
            return True
        except Exception as e:
            logger.error(f"snpEff annotation failed: {e}")
            return False

    def annotate_with_clinvar(self, variants: List[Dict]) -> List[Dict]:
        """ClinVar annotation 추가"""
        if not self.clinvar_vcf or not os.path.exists(self.clinvar_vcf):
            logger.warning("ClinVar VCF not found, skipping ClinVar annotation")
            return variants

        try:
            import pysam
            clinvar_tbx = pysam.VariantFile(self.clinvar_vcf)

            for var in variants:
                chrom = var.get("chrom", "")
                pos = var.get("pos", 0)
                ref = var.get("ref", "")
                alt = var.get("alt", "")

                try:
                    for rec in clinvar_tbx.fetch(chrom, max(0, pos - 1), pos + len(ref)):
                        if rec.pos == pos and rec.ref == ref:
                            for alt_allele in rec.alts or []:
                                if alt_allele == alt:
                                    clnsig = rec.info.get("CLNSIG", None)
                                    clndn = rec.info.get("CLNDN", None)
                                    clnrevstat = rec.info.get("CLNREVSTAT", None)
                                    var["clinvar"] = {
                                        "significance": _join_tuple(clnsig),
                                        "disease": _join_tuple(clndn),
                                        "review_status": _join_tuple(clnrevstat),
                                        "id": rec.id
                                    }
                                    break
                except Exception:
                    pass

            clinvar_tbx.close()
            logger.info(f"ClinVar annotation complete for {len(variants)} variants")

        except ImportError:
            logger.warning("pysam not available, skipping ClinVar annotation")
        except Exception as e:
            logger.error(f"ClinVar annotation error: {e}")

        return variants

    def annotate_with_gnomad(self, variants: List[Dict]) -> List[Dict]:
        """gnomAD allele frequency annotation 추가"""
        if not self.gnomad_vcf or not os.path.exists(self.gnomad_vcf):
            logger.warning("gnomAD VCF not found, skipping gnomAD annotation")
            return variants

        try:
            import pysam
            gnomad_tbx = pysam.VariantFile(self.gnomad_vcf)

            for var in variants:
                chrom = var.get("chrom", "")
                pos = var.get("pos", 0)
                ref = var.get("ref", "")
                alt = var.get("alt", "")

                try:
                    for rec in gnomad_tbx.fetch(chrom, max(0, pos - 1), pos + len(ref)):
                        if rec.pos == pos and rec.ref == ref:
                            for i, alt_allele in enumerate(rec.alts or []):
                                if alt_allele == alt:
                                    af = rec.info.get("AF", None)
                                    if af and isinstance(af, tuple) and len(af) > i:
                                        af_val = af[i]
                                    elif af:
                                        af_val = af[0] if isinstance(af, tuple) else af
                                    else:
                                        af_val = None

                                    var["gnomad"] = {
                                        "af": af_val,
                                        "af_popmax": rec.info.get("AF_popmax", None),
                                        "nhomalt": rec.info.get("nhomalt", None)
                                    }
                                    break
                except Exception:
                    pass

            gnomad_tbx.close()
            logger.info(f"gnomAD annotation complete for {len(variants)} variants")

        except ImportError:
            logger.warning("pysam not available, skipping gnomAD annotation")
        except Exception as e:
            logger.error(f"gnomAD annotation error: {e}")

        return variants

    def annotate_with_dbsnp(self, variants: List[Dict]) -> List[Dict]:
        """dbSNP rsID annotation 추가"""
        if not self.dbsnp_vcf or not os.path.exists(self.dbsnp_vcf):
            logger.warning("dbSNP VCF not found, skipping dbSNP annotation")
            return variants

        try:
            import pysam
            dbsnp_tbx = pysam.VariantFile(self.dbsnp_vcf)

            for var in variants:
                chrom = var.get("chrom", "")
                pos = var.get("pos", 0)

                try:
                    for rec in dbsnp_tbx.fetch(chrom, max(0, pos - 1), pos + 1):
                        if rec.pos == pos:
                            var["dbsnp"] = {"rsid": rec.id}
                            break
                except Exception:
                    pass

            dbsnp_tbx.close()
            logger.info(f"dbSNP annotation complete for {len(variants)} variants")

        except ImportError:
            logger.warning("pysam not available, skipping dbSNP annotation")
        except Exception as e:
            logger.error(f"dbSNP annotation error: {e}")

        return variants


def _join_tuple(val) -> Optional[str]:
    """pysam tuple 값을 문자열로 변환"""
    if val is None:
        return None
    if isinstance(val, (tuple, list)):
        return ",".join(str(v) for v in val)
    return str(val)


# ─── Review Page Generator ─────────────────────────────────

class ReviewPageGenerator:
    """
    리뷰 페이지용 JSON 및 관련 파일을 생성합니다.
    
    출력 구조:
    output_dir/
    ├── review_data.json          # 리뷰 페이지 메인 데이터
    ├── variants_annotated.json   # Annotation된 변이 목록
    ├── summary/                  # 요약 리포트
    ├── snapshots/                # IGV 스냅샷 이미지
    └── metadata.json             # 메타데이터
    """

    def generate_review_data(
        self,
        job: Job,
        summary_data: Dict[str, Any],
        annotated_variants: List[Dict],
        snapshot_files: List[str]
    ) -> Dict[str, Any]:
        """
        리뷰 페이지용 JSON 데이터를 생성합니다.
        """
        review_data = {
            "version": "1.0",
            "service_code": "carrier_screening",
            "order_id": job.order_id,
            "sample_name": job.sample_name,
            "generated_at": _now_iso(),
            
            # 요약 정보
            "summary": {
                "total_variants": len(annotated_variants),
                "pathogenic_count": sum(
                    1 for v in annotated_variants
                    if _is_pathogenic(v)
                ),
                "likely_pathogenic_count": sum(
                    1 for v in annotated_variants
                    if _is_likely_pathogenic(v)
                ),
                "vus_count": sum(
                    1 for v in annotated_variants
                    if _is_vus(v)
                ),
                "pipeline_summary": summary_data
            },

            # 변이 목록 (annotation 포함)
            "variants": annotated_variants,

            # 스냅샷 이미지 참조
            "snapshots": [
                {
                    "filename": os.path.basename(f),
                    "path": f"snapshots/{os.path.basename(f)}"
                }
                for f in snapshot_files
            ],

            # CNV 결과 (있는 경우)
            "cnv_results": summary_data.get("cnv", {}),

            # SV 결과 (있는 경우)
            "sv_results": summary_data.get("sv", {}),

            # Repeat Expansion 결과
            "repeat_expansion": summary_data.get("repeat_expansion", {}),

            # Pseudogene 결과
            "pseudogene": summary_data.get("pseudogene", {}),

            # 리뷰 상태 (Platform에서 관리)
            "review_status": "PENDING",
            "reviewer": None,
            "review_notes": None
        }

        return review_data

    def generate_metadata(self, job: Job) -> Dict[str, Any]:
        """메타데이터 JSON 생성"""
        return {
            "order_id": job.order_id,
            "service_code": "carrier_screening",
            "sample_name": job.sample_name,
            "work_dir": job.work_dir,
            "pipeline_version": _get_pipeline_version(),
            "annotation_databases": {
                "clinvar": settings.clinvar_vcf or "not_configured",
                "gnomad": settings.gnomad_vcf or "not_configured",
                "dbsnp": settings.dbsnp_vcf or "not_configured",
                "snpeff": settings.snpeff_db
            },
            "generated_at": _now_iso()
        }


# ─── VCF Parser ────────────────────────────────────────────

def parse_vcf_variants(vcf_path: str) -> List[Dict]:
    """VCF 파일에서 변이를 파싱합니다."""
    variants = []

    opener = gzip.open if vcf_path.endswith(".gz") else open
    mode = "rt" if vcf_path.endswith(".gz") else "r"

    try:
        with opener(vcf_path, mode) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue

                chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]

                # INFO 필드 파싱
                info_dict = {}
                for item in info.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info_dict[k] = v
                    else:
                        info_dict[item] = True

                # 각 ALT allele에 대해 변이 생성
                for alt_allele in alt.split(","):
                    variant = {
                        "chrom": chrom,
                        "pos": int(pos),
                        "id": vid if vid != "." else None,
                        "ref": ref,
                        "alt": alt_allele,
                        "qual": float(qual) if qual != "." else None,
                        "filter": filt,
                        "info": info_dict,
                        # snpEff annotation (ANN 필드)
                        "gene": _extract_gene_from_info(info_dict),
                        "effect": _extract_effect_from_info(info_dict),
                        "hgvsc": _extract_hgvsc_from_info(info_dict),
                        "hgvsp": _extract_hgvsp_from_info(info_dict),
                        # Genotype (있는 경우)
                        "genotype": parts[9].split(":")[0] if len(parts) > 9 else None
                    }
                    variants.append(variant)

        logger.info(f"Parsed {len(variants)} variants from {vcf_path}")

    except Exception as e:
        logger.error(f"Failed to parse VCF {vcf_path}: {e}")

    return variants


def _extract_gene_from_info(info: Dict) -> Optional[str]:
    """INFO 필드에서 유전자명 추출"""
    # snpEff ANN 필드
    ann = info.get("ANN", "")
    if ann and isinstance(ann, str):
        parts = ann.split("|")
        if len(parts) > 3:
            return parts[3]
    return info.get("GENE", info.get("Gene", None))


def _extract_effect_from_info(info: Dict) -> Optional[str]:
    """INFO 필드에서 변이 효과 추출"""
    ann = info.get("ANN", "")
    if ann and isinstance(ann, str):
        parts = ann.split("|")
        if len(parts) > 1:
            return parts[1]
    return None


def _extract_hgvsc_from_info(info: Dict) -> Optional[str]:
    """INFO 필드에서 HGVS.c 추출"""
    ann = info.get("ANN", "")
    if ann and isinstance(ann, str):
        parts = ann.split("|")
        if len(parts) > 9:
            return parts[9] if parts[9] else None
    return None


def _extract_hgvsp_from_info(info: Dict) -> Optional[str]:
    """INFO 필드에서 HGVS.p 추출"""
    ann = info.get("ANN", "")
    if ann and isinstance(ann, str):
        parts = ann.split("|")
        if len(parts) > 10:
            return parts[10] if parts[10] else None
    return None


def _is_pathogenic(variant: Dict) -> bool:
    """ClinVar Pathogenic 여부"""
    clinvar = variant.get("clinvar", {})
    sig = (clinvar.get("significance") or "").lower()
    return "pathogenic" in sig and "likely" not in sig


def _is_likely_pathogenic(variant: Dict) -> bool:
    """ClinVar Likely Pathogenic 여부"""
    clinvar = variant.get("clinvar", {})
    sig = (clinvar.get("significance") or "").lower()
    return "likely_pathogenic" in sig or "likely pathogenic" in sig


def _is_vus(variant: Dict) -> bool:
    """VUS (Variant of Uncertain Significance) 여부"""
    clinvar = variant.get("clinvar", {})
    sig = (clinvar.get("significance") or "").lower()
    return "uncertain" in sig


def _now_iso() -> str:
    from datetime import datetime
    return datetime.now().isoformat()


def _get_pipeline_version() -> str:
    """파이프라인 버전 조회"""
    try:
        version_file = os.path.join(
            settings.carrier_screening_pipeline_dir, "VERSION"
        )
        if os.path.exists(version_file):
            with open(version_file) as f:
                return f.read().strip()
    except Exception:
        pass
    return "unknown"


# ─── Carrier Screening Plugin ──────────────────────────────

class CarrierScreeningPlugin(ServicePlugin):
    """Carrier Screening 서비스 플러그인"""

    def __init__(self):
        self._annotator = VCFAnnotator()
        self._review_generator = ReviewPageGenerator()

    @property
    def service_code(self) -> str:
        return "carrier_screening"

    @property
    def display_name(self) -> str:
        return "Carrier Screening"

    async def prepare_inputs(self, job: Job) -> bool:
        """FASTQ 파일 준비 및 디렉토리 구조 생성"""
        try:
            # 디렉토리 생성
            for dir_path in [job.fastq_dir, job.analysis_dir, job.output_dir, job.log_dir]:
                if dir_path:
                    os.makedirs(dir_path, exist_ok=True)

            # URL이 제공된 경우 다운로드
            if job.fastq_r1_url and job.fastq_r2_url:
                for url, label in [
                    (job.fastq_r1_url, "R1"), (job.fastq_r2_url, "R2")
                ]:
                    dest = os.path.join(job.fastq_dir, os.path.basename(url))
                    if not os.path.exists(dest):
                        logger.info(
                            f"[carrier_screening] Downloading {label}: {url} -> {dest}"
                        )
                        subprocess.run(
                            ["wget", "-q", "-O", dest, url],
                            check=True, timeout=7200
                        )
                    else:
                        logger.info(f"[carrier_screening] {label} already exists: {dest}")

            # 로컬 경로가 제공된 경우 심볼릭 링크
            elif job.fastq_r1_path and job.fastq_r2_path:
                for src in [job.fastq_r1_path, job.fastq_r2_path]:
                    if not os.path.exists(src):
                        logger.error(f"[carrier_screening] FASTQ not found: {src}")
                        return False
                    dst = os.path.join(job.fastq_dir, os.path.basename(src))
                    if not os.path.exists(dst):
                        os.symlink(os.path.abspath(src), dst)
                        logger.info(f"[carrier_screening] Linked {src} -> {dst}")

            # FASTQ 파일 존재 확인
            fastq_files = (
                glob.glob(os.path.join(job.fastq_dir, "*.fastq.gz")) +
                glob.glob(os.path.join(job.fastq_dir, "*.fq.gz"))
            )
            if not fastq_files:
                logger.error(
                    f"[carrier_screening] No FASTQ files found in {job.fastq_dir}"
                )
                return False

            logger.info(
                f"[carrier_screening] Input ready: {len(fastq_files)} FASTQ files "
                f"in {job.fastq_dir}"
            )
            return True

        except Exception as e:
            logger.error(f"[carrier_screening] Failed to prepare inputs: {e}")
            return False

    async def get_pipeline_command(self, job: Job) -> str:
        """Nextflow 파이프라인 실행 명령어 생성"""
        pipeline_dir = settings.carrier_screening_pipeline_dir
        nf_script = os.path.join(pipeline_dir, "bin", "main.nf")
        analysis_dir = job.analysis_dir
        output_dir = job.output_dir

        # 기본 명령어
        cmd_parts = [
            settings.nextflow_executable,
            "run", nf_script,
            f"--fastq_dir {job.fastq_dir}",
            f"--outdir {analysis_dir}",
            f"--output_dir {output_dir}",
            f"--sample_name {job.sample_name}",
            f"-work-dir {analysis_dir}/work",
            "-with-trace",
            "-with-timeline",
            "-with-report"
        ]

        # Nextflow 설정 파일
        nf_config = settings.nextflow_config or os.path.join(pipeline_dir, "nextflow.config")
        if os.path.exists(nf_config):
            cmd_parts.append(f"-c {nf_config}")

        # Reference genome
        if settings.ref_fasta:
            cmd_parts.append(f"--ref_fasta {settings.ref_fasta}")
        if settings.ref_fai:
            cmd_parts.append(f"--ref_fai {settings.ref_fai}")
        if settings.ref_dict:
            cmd_parts.append(f"--ref_dict {settings.ref_dict}")
        if settings.ref_bwa_indices:
            cmd_parts.append(f"--ref_bwa_indices {settings.ref_bwa_indices}")

        # 서비스별 추가 파라미터 (job.params에서 동적으로)
        for key, value in job.params.items():
            if value is not None and value != "":
                cmd_parts.append(f"--{key} {value}")

        command = " ".join(cmd_parts)
        logger.info(f"[carrier_screening] Pipeline command: {command}")
        return command

    async def check_completion(self, job: Job) -> bool:
        """파이프라인 완료 확인"""
        analysis_dir = job.analysis_dir
        output_dir = job.output_dir

        # 1. Summary report 확인
        summary_files = (
            glob.glob(os.path.join(analysis_dir, "**", "*summary*report*"), recursive=True) +
            glob.glob(os.path.join(output_dir, "**", "*summary*"), recursive=True)
        )

        # 2. VCF 결과 확인
        vcf_files = (
            glob.glob(os.path.join(analysis_dir, "**", "*.vcf.gz"), recursive=True) +
            glob.glob(os.path.join(output_dir, "**", "*.vcf.gz"), recursive=True)
        )

        if summary_files or vcf_files:
            logger.info(
                f"[carrier_screening] Completion check passed: "
                f"{len(summary_files)} summary files, {len(vcf_files)} VCF files"
            )
            return True

        # 3. Nextflow trace 파일에서 성공 확인
        trace_file = os.path.join(analysis_dir, "work", "trace.txt")
        if os.path.exists(trace_file):
            try:
                with open(trace_file) as f:
                    content = f.read()
                    if "COMPLETED" in content:
                        logger.info(
                            "[carrier_screening] Completion confirmed via trace file"
                        )
                        return True
            except Exception:
                pass

        logger.warning(
            f"[carrier_screening] Completion check failed for {job.order_id}"
        )
        return False

    async def process_results(self, job: Job) -> bool:
        """
        결과 후처리: Annotation + 리뷰 페이지 JSON 생성.
        
        1. VCF에서 변이 파싱
        2. ClinVar, gnomAD, dbSNP annotation
        3. 리뷰 페이지 JSON 생성
        4. 스냅샷 이미지 정리
        5. 메타데이터 생성
        """
        try:
            output_dir = job.output_dir
            analysis_dir = job.analysis_dir
            os.makedirs(output_dir, exist_ok=True)

            # ── Step 1: VCF 파일 찾기 및 변이 파싱 ──
            vcf_files = self._find_vcf_files(analysis_dir, output_dir)
            all_variants = []
            for vcf_path in vcf_files:
                variants = parse_vcf_variants(vcf_path)
                all_variants.extend(variants)

            logger.info(
                f"[carrier_screening] Parsed {len(all_variants)} total variants "
                f"from {len(vcf_files)} VCF files"
            )

            # ── Step 2: Annotation ──
            annotated_variants = self._annotator.annotate_with_clinvar(all_variants)
            annotated_variants = self._annotator.annotate_with_gnomad(annotated_variants)
            annotated_variants = self._annotator.annotate_with_dbsnp(annotated_variants)

            # Annotation된 변이 저장
            annotated_path = os.path.join(output_dir, "variants_annotated.json")
            with open(annotated_path, "w", encoding="utf-8") as f:
                json.dump(annotated_variants, f, indent=2, ensure_ascii=False, default=str)
            logger.info(f"[carrier_screening] Saved annotated variants: {annotated_path}")

            # ── Step 3: Summary 데이터 수집 ──
            summary_data = self._collect_summary_data(analysis_dir, output_dir)

            # ── Step 4: 스냅샷 이미지 정리 ──
            snapshot_dir = os.path.join(output_dir, "snapshots")
            os.makedirs(snapshot_dir, exist_ok=True)
            snapshot_files = self._collect_snapshots(analysis_dir, output_dir, snapshot_dir)

            # ── Step 5: 리뷰 페이지 JSON 생성 ──
            review_data = self._review_generator.generate_review_data(
                job, summary_data, annotated_variants, snapshot_files
            )
            review_path = os.path.join(output_dir, "review_data.json")
            with open(review_path, "w", encoding="utf-8") as f:
                json.dump(review_data, f, indent=2, ensure_ascii=False, default=str)
            logger.info(f"[carrier_screening] Saved review data: {review_path}")

            # ── Step 6: 메타데이터 생성 ──
            metadata = self._review_generator.generate_metadata(job)
            metadata_path = os.path.join(output_dir, "metadata.json")
            with open(metadata_path, "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=2, ensure_ascii=False)
            logger.info(f"[carrier_screening] Saved metadata: {metadata_path}")

            # ── Step 7: Summary 리포트 복사 ──
            summary_dir = os.path.join(output_dir, "summary")
            os.makedirs(summary_dir, exist_ok=True)
            for pattern in ["*summary*report*", "*detailed*report*"]:
                for f in glob.glob(
                    os.path.join(analysis_dir, "**", pattern), recursive=True
                ):
                    dst = os.path.join(summary_dir, os.path.basename(f))
                    if not os.path.exists(dst):
                        shutil.copy2(f, dst)

            logger.info(
                f"[carrier_screening] Result processing complete for {job.order_id}"
            )
            return True

        except Exception as e:
            logger.error(
                f"[carrier_screening] Result processing failed for {job.order_id}: {e}",
                exc_info=True
            )
            return False

    async def get_output_files(self, job: Job) -> List[OutputFile]:
        """플랫폼에 업로드할 출력 파일 목록"""
        output_dir = job.output_dir
        files = []

        # 리뷰 데이터 JSON (최우선)
        review_json = os.path.join(output_dir, "review_data.json")
        if os.path.exists(review_json):
            files.append(OutputFile(
                file_path=review_json,
                file_type="review_json",
                file_name="review_data.json",
                content_type="application/json"
            ))

        # Annotation된 변이 JSON
        variants_json = os.path.join(output_dir, "variants_annotated.json")
        if os.path.exists(variants_json):
            files.append(OutputFile(
                file_path=variants_json,
                file_type="variants_json",
                file_name="variants_annotated.json",
                content_type="application/json"
            ))

        # 메타데이터
        metadata_json = os.path.join(output_dir, "metadata.json")
        if os.path.exists(metadata_json):
            files.append(OutputFile(
                file_path=metadata_json,
                file_type="metadata",
                file_name="metadata.json",
                content_type="application/json"
            ))

        # 스냅샷 이미지
        snapshot_dir = os.path.join(output_dir, "snapshots")
        if os.path.exists(snapshot_dir):
            for img_file in sorted(glob.glob(os.path.join(snapshot_dir, "*.png"))):
                files.append(OutputFile(
                    file_path=img_file,
                    file_type="snapshot",
                    content_type="image/png"
                ))
            for img_file in sorted(glob.glob(os.path.join(snapshot_dir, "*.svg"))):
                files.append(OutputFile(
                    file_path=img_file,
                    file_type="snapshot",
                    content_type="image/svg+xml"
                ))

        # Summary 리포트
        summary_dir = os.path.join(output_dir, "summary")
        if os.path.exists(summary_dir):
            for f in sorted(glob.glob(os.path.join(summary_dir, "*"))):
                files.append(OutputFile(
                    file_path=f,
                    file_type="summary_report",
                    content_type="application/octet-stream"
                ))

        # VCF 파일
        for vcf_file in glob.glob(os.path.join(output_dir, "**", "*.vcf.gz"), recursive=True):
            files.append(OutputFile(
                file_path=vcf_file,
                file_type="vcf",
                content_type="application/gzip"
            ))

        logger.info(
            f"[carrier_screening] Output files: {len(files)} files for upload"
        )
        return files

    def get_progress_stages(self) -> Dict[str, int]:
        """Carrier Screening 파이프라인 진행률 단계"""
        return {
            "ALIGN_AND_SORT": 15,
            "MARK_DUPLICATES": 20,
            "CALL_VARIANTS": 35,
            "COLLECT_READ_COUNTS": 40,
            "GCNV": 50,
            "POSTPROCESS_GCNV": 55,
            "MANTA_SV": 60,
            "EXPANSION_HUNTER": 65,
            "PARAPHASE_RUN": 70,
            "SMACA_RUN": 72,
            "DEPTH_ANALYSIS": 75,
            "FALLBACK_ANALYSIS": 78,
            "GENERATE_VISUAL_EVIDENCE": 82,
            "GENERATE_SUMMARY_REPORT": 85
        }

    def validate_params(self, params: Dict[str, Any]) -> tuple:
        """Carrier Screening 파라미터 유효성 검사"""
        # backbone_bed는 권장 파라미터
        if "backbone_bed" in params:
            bed_path = params["backbone_bed"]
            if not os.path.exists(bed_path):
                return False, f"backbone_bed not found: {bed_path}"
        return True, ""

    # ─── Private Helpers ───────────────────────────────────

    def _find_vcf_files(self, analysis_dir: str, output_dir: str) -> List[str]:
        """분석 결과에서 VCF 파일을 찾습니다."""
        vcf_files = []

        # 우선순위: filtered VCF > 일반 VCF
        search_dirs = [output_dir, analysis_dir]
        patterns = [
            "*_filtered*.vcf.gz",
            "*_filtered*.vcf",
            "*.vcf.gz",
            "*.vcf"
        ]

        found_set = set()
        for search_dir in search_dirs:
            for pattern in patterns:
                for f in glob.glob(
                    os.path.join(search_dir, "**", pattern), recursive=True
                ):
                    if f not in found_set:
                        vcf_files.append(f)
                        found_set.add(f)

        logger.info(f"[carrier_screening] Found {len(vcf_files)} VCF files")
        return vcf_files

    def _collect_summary_data(
        self, analysis_dir: str, output_dir: str
    ) -> Dict[str, Any]:
        """분석 결과에서 요약 데이터를 수집합니다."""
        summary = {}

        # Summary report JSON/TXT 파일 찾기
        for search_dir in [output_dir, analysis_dir]:
            for pattern in ["*summary*report*.json", "*summary*report*.txt"]:
                for f in glob.glob(
                    os.path.join(search_dir, "**", pattern), recursive=True
                ):
                    try:
                        if f.endswith(".json"):
                            with open(f) as fh:
                                data = json.load(fh)
                                if isinstance(data, dict):
                                    summary.update(data)
                        else:
                            with open(f) as fh:
                                summary["summary_text"] = fh.read()
                    except Exception as e:
                        logger.warning(f"Failed to read summary file {f}: {e}")

        # Paraphase 결과
        for f in glob.glob(
            os.path.join(analysis_dir, "**", "*paraphase*.json"), recursive=True
        ):
            try:
                with open(f) as fh:
                    summary["pseudogene"] = json.load(fh)
                break
            except Exception:
                pass

        # Expansion Hunter 결과
        for f in glob.glob(
            os.path.join(analysis_dir, "**", "*eh*.json"), recursive=True
        ):
            try:
                with open(f) as fh:
                    summary["repeat_expansion"] = json.load(fh)
                break
            except Exception:
                pass

        return summary

    def _collect_snapshots(
        self, analysis_dir: str, output_dir: str, snapshot_dir: str
    ) -> List[str]:
        """스냅샷 이미지를 수집하여 output 디렉토리로 복사합니다."""
        snapshot_files = []

        for search_dir in [output_dir, analysis_dir]:
            for pattern in ["*.png", "*.svg", "*.jpg"]:
                for f in glob.glob(
                    os.path.join(search_dir, "**", "snapshots", pattern),
                    recursive=True
                ):
                    dst = os.path.join(snapshot_dir, os.path.basename(f))
                    if not os.path.exists(dst):
                        shutil.copy2(f, dst)
                    snapshot_files.append(dst)

                # 시각화 결과 (IGV 등)
                for f in glob.glob(
                    os.path.join(search_dir, "**", "visual*", pattern),
                    recursive=True
                ):
                    dst = os.path.join(snapshot_dir, os.path.basename(f))
                    if not os.path.exists(dst):
                        shutil.copy2(f, dst)
                    snapshot_files.append(dst)

        logger.info(
            f"[carrier_screening] Collected {len(snapshot_files)} snapshot files"
        )
        return snapshot_files


def create_plugin() -> ServicePlugin:
    """플러그인 팩토리 함수"""
    return CarrierScreeningPlugin()
