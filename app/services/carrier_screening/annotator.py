"""
Variant Annotator Module

VCF 변이에 대해 ClinVar, gnomAD, ClinGen, MANE, dbSNP 등의 annotation을 수행합니다.
phenotype_portal/main.py의 annotation 로직을 모듈화한 것입니다.

주요 기능:
    - ClinVar 로컬 VCF 조회 (significance, review status, conflict 분석)
    - gnomAD 로컬 VCF 조회 (exomes + genomes AF)
    - ClinGen 유전자 용량 민감도 조회
    - MANE RefSeq 표준 전사체 매핑
    - dbSNP URL 생성
"""

import os
import re
import csv
import gzip
import logging
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List
from collections import Counter, defaultdict

logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════
# ClinVar Lookup
# ══════════════════════════════════════════════════════════════

class ClinVarAnnotator:
    """로컬 ClinVar VCF를 사용한 변이 annotation"""

    def __init__(self, clinvar_vcf_path: str):
        self._vcf_path = clinvar_vcf_path
        self._has_chr: Optional[bool] = None

        if clinvar_vcf_path and os.path.exists(clinvar_vcf_path):
            logger.info(f"ClinVar VCF: {clinvar_vcf_path}")
        else:
            logger.warning(f"ClinVar VCF not found: {clinvar_vcf_path}")

    def _detect_has_chr(self) -> Optional[bool]:
        try:
            import pysam
            vf = pysam.VariantFile(self._vcf_path)
            contigs = list(vf.header.contigs)
            vf.close()
            if not contigs:
                return None
            return any(str(c).startswith("chr") for c in contigs)
        except Exception:
            return None

    def _normalize_chrom(self, chrom: str) -> str:
        if self._has_chr is None:
            self._has_chr = self._detect_has_chr()
            if self._has_chr is None:
                self._has_chr = True

        chrom = str(chrom)
        if self._has_chr:
            return chrom if chrom.startswith("chr") else "chr" + chrom
        return chrom[3:] if chrom.startswith("chr") else chrom

    @staticmethod
    def _first_val(info, key: str) -> str:
        if key not in info:
            return ""
        v = info[key]
        if v is None:
            return ""
        if isinstance(v, (list, tuple)):
            return str(v[0]) if v else ""
        return str(v)

    @staticmethod
    def _all_vals_joined(info, key: str, sep: str = "|") -> str:
        if key not in info:
            return ""
        v = info[key]
        if v is None:
            return ""
        if isinstance(v, (list, tuple)):
            vals = [str(x) for x in v if x not in (None, "", ".")]
            return sep.join(vals)
        s = str(v).strip()
        return "" if s in ("", ".") else s

    @staticmethod
    def _sig_primary(clnsig_raw: str) -> str:
        if not clnsig_raw:
            return ""
        s = clnsig_raw.replace("_", " ").strip().lower()
        if "conflict" in s:
            return "Conflicting"
        if "pathogenic" in s and "likely" in s:
            return "Likely pathogenic"
        if "pathogenic" in s:
            return "Pathogenic"
        if "uncertain" in s or "vus" in s:
            return "VUS"
        if "benign" in s and "likely" in s:
            return "Likely benign"
        if "benign" in s:
            return "Benign"
        return clnsig_raw.strip()

    @staticmethod
    def _is_conflicting(clnsig_raw: str) -> bool:
        if not clnsig_raw:
            return False
        s = clnsig_raw.lower()
        return ("conflicting" in s) or ("conflict" in s) or (("benign" in s) and ("pathogenic" in s))

    @staticmethod
    def _review_stars(revstat_raw: str) -> int:
        if not revstat_raw:
            return 0
        s = revstat_raw.replace("_", " ").lower()
        if "practice guideline" in s:
            return 4
        if "reviewed by expert panel" in s:
            return 3
        if "multiple submitters" in s and "no conflicts" in s:
            return 2
        if "single submitter" in s:
            return 1
        return 0

    @staticmethod
    def _normalize_sig_label(label: str) -> str:
        l = (label or "").strip().replace("_", " ")
        low = l.lower()
        if "likely benign" in low:
            return "Likely benign"
        if low == "benign":
            return "Benign"
        if "likely pathogenic" in low:
            return "Likely pathogenic"
        if low == "pathogenic":
            return "Pathogenic"
        if "uncertain" in low or "vus" in low:
            return "VUS"
        if "conflict" in low:
            return "Conflicting"
        return l.strip()

    @staticmethod
    def _parse_clnsigconf_counts(clnsigconf_raw: str) -> Counter:
        c = Counter()
        if not clnsigconf_raw:
            return c
        parts = re.split(r"[|,]", clnsigconf_raw)
        for p in parts:
            p = p.strip()
            if not p:
                continue
            m = re.match(r"^(.*?)(?:\s*\(\s*(\d+)\s*\))?$", p)
            if not m:
                continue
            label = ClinVarAnnotator._normalize_sig_label(m.group(1))
            n = int(m.group(2)) if m.group(2) else 1
            if label:
                c[label] += n
        if "Conflicting" in c:
            del c["Conflicting"]
        return c

    @staticmethod
    def _format_conflict_counts(counter: Counter) -> str:
        if not counter:
            return ""
        benignish = ["Benign", "Likely benign"]
        pathish = ["Pathogenic", "Likely pathogenic"]

        def fmt_group(keys):
            items = [f"{k} ({counter[k]})" for k in keys if counter.get(k, 0) > 0]
            return " / ".join(items)

        left = fmt_group(benignish)
        right = fmt_group(pathish)
        if left and right:
            return f"{left} vs {right}"
        return " / ".join(f"{k} ({v})" for k, v in counter.most_common())

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[Dict[str, Any]]:
        """ClinVar에서 변이를 조회합니다."""
        if not self._vcf_path or not os.path.exists(self._vcf_path):
            return None

        try:
            import pysam
            cv_chrom = self._normalize_chrom(chrom)
            cv = pysam.VariantFile(self._vcf_path)

            for rec in cv.fetch(cv_chrom, pos - 1, pos):
                if rec.pos != pos:
                    continue
                if rec.ref != ref:
                    continue
                if not rec.alts or alt not in rec.alts:
                    continue

                info = rec.info
                clnsig = self._all_vals_joined(info, "CLNSIG", sep="|")
                clnsigconf = self._all_vals_joined(info, "CLNSIGCONF", sep="|")
                revstat = self._first_val(info, "CLNREVSTAT")
                clndn = self._first_val(info, "CLNDN")
                variation_id = rec.id or ""

                rs_raw = self._first_val(info, "RS")
                rsid = ""
                if rs_raw and rs_raw not in (".", "0"):
                    m = re.search(r"(\d+)", str(rs_raw))
                    if m:
                        rsid = f"rs{m.group(1)}"

                conflict_detail = ""
                if self._is_conflicting(clnsig):
                    counts = self._parse_clnsigconf_counts(clnsigconf)
                    conflict_detail = self._format_conflict_counts(counts)

                result = {
                    "clnsig": clnsig,
                    "clnsig_primary": self._sig_primary(clnsig),
                    "conflicting": self._is_conflicting(clnsig),
                    "conflict_summary": conflict_detail,
                    "revstat": revstat,
                    "stars": self._review_stars(revstat),
                    "clndn": clndn,
                    "clnvc": self._first_val(info, "CLNVC"),
                    "allele_id": self._first_val(info, "ALLELEID"),
                    "variation_id": variation_id,
                    "rsid": rsid,
                }
                cv.close()
                return result

            cv.close()
            return None

        except Exception as e:
            logger.error(f"ClinVar lookup error {chrom}:{pos} {ref}>{alt}: {e}")
            return None


# ══════════════════════════════════════════════════════════════
# gnomAD Lookup
# ══════════════════════════════════════════════════════════════

class GnomADAnnotator:
    """로컬 gnomAD VCF를 사용한 allele frequency 조회"""

    _AF_KEYS = ["AF", "AF_total", "AF_joint", "AF_popmax", "AF_POPMAX"]

    def __init__(self, gnomad_dir: str,
                 genomes_glob: str = "gnomad.genomes.v*.sites*.bgz",
                 exomes_glob: str = "gnomad.exomes.v*.sites*.bgz"):
        self._gnomad_dir = gnomad_dir
        self._genomes_glob = genomes_glob
        self._exomes_glob = exomes_glob
        self._vcf_cache: Dict[str, Any] = {}
        self._has_chr_cache: Dict[str, bool] = {}
        self._pick_cache: Dict[Tuple[str, str], Optional[str]] = {}

        if gnomad_dir and os.path.isdir(gnomad_dir):
            logger.info(f"gnomAD directory: {gnomad_dir}")
        else:
            logger.warning(f"gnomAD directory not found: {gnomad_dir}")

    @staticmethod
    def _first_numeric(val) -> Optional[float]:
        if val is None:
            return None
        if isinstance(val, (list, tuple)):
            val = val[0] if val else None
        if val is None:
            return None
        try:
            return float(val)
        except Exception:
            try:
                return float(str(val))
            except Exception:
                return None

    def _detect_vcf_has_chr(self, vcf_path: str) -> bool:
        try:
            import pysam
            vf = pysam.VariantFile(vcf_path)
            contigs = list(vf.header.contigs)
            vf.close()
            return any(str(c).startswith("chr") for c in contigs)
        except Exception:
            return True

    def _norm_chrom(self, chrom: str, vcf_has_chr: bool) -> str:
        chrom = str(chrom)
        if vcf_has_chr:
            return chrom if chrom.startswith("chr") else "chr" + chrom
        return chrom[3:] if chrom.startswith("chr") else chrom

    def _open_vcf(self, vcf_path: str):
        import pysam
        vf = self._vcf_cache.get(vcf_path)
        if vf is None:
            vf = pysam.VariantFile(vcf_path)
            self._vcf_cache[vcf_path] = vf
            self._has_chr_cache[vcf_path] = self._detect_vcf_has_chr(vcf_path)
        return vf

    def _pick_file(self, dataset: str, chrom: str) -> Optional[str]:
        key = (dataset, str(chrom))
        if key in self._pick_cache:
            return self._pick_cache[key]

        if not self._gnomad_dir or not os.path.isdir(self._gnomad_dir):
            self._pick_cache[key] = None
            return None

        c = str(chrom)
        c_nochr = c[3:] if c.startswith("chr") else c
        token_chr = f"chr{c_nochr}".lower()

        glob_pat = self._exomes_glob if dataset == "exomes" else self._genomes_glob
        candidates = list(Path(self._gnomad_dir).glob(glob_pat))

        best = None
        for p in candidates:
            name = p.name.lower()
            if token_chr in name:
                best = str(p)
                break

        if best is None:
            token_no = c_nochr.lower()
            for p in candidates:
                name = p.name.lower()
                if f".{token_no}." in name or f"_{token_no}_" in name:
                    best = str(p)
                    break

        if best is None and len(candidates) == 1:
            best = str(candidates[0])

        self._pick_cache[key] = best
        return best

    def _extract_alt_af(self, rec, alt: str, key: str) -> Optional[float]:
        if key not in rec.info:
            return None
        val = rec.info[key]
        if isinstance(val, (list, tuple)):
            if not rec.alts:
                return self._first_numeric(val)
            try:
                idx = list(rec.alts).index(alt)
            except ValueError:
                idx = 0
            if idx < len(val):
                return self._first_numeric(val[idx])
            return self._first_numeric(val[0])
        return self._first_numeric(val)

    def _lookup_one(self, dataset: str, chrom: str, pos: int, ref: str, alt: str) -> Optional[float]:
        vcf_path = self._pick_file(dataset, chrom)
        if not vcf_path or not os.path.exists(vcf_path):
            return None

        has_index = (os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi"))
        if not has_index:
            return None

        try:
            vf = self._open_vcf(vcf_path)
            vcf_has_chr = self._has_chr_cache.get(vcf_path, True)
            qchrom = self._norm_chrom(chrom, vcf_has_chr)

            for r in vf.fetch(qchrom, pos - 1, pos):
                if r.pos != pos:
                    continue
                if r.ref != ref:
                    continue
                if not r.alts or alt not in r.alts:
                    continue

                for k in self._AF_KEYS:
                    af = self._extract_alt_af(r, alt, k)
                    if af is not None:
                        return af

        except Exception as e:
            logger.debug(f"gnomAD lookup error {dataset} {chrom}:{pos}: {e}")
            return None

        return None

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
        """gnomAD에서 allele frequency를 조회합니다."""
        out = {"af": None, "exomes_af": None, "genomes_af": None, "source": ""}

        if not self._gnomad_dir or not os.path.isdir(self._gnomad_dir):
            return out

        ex_af = self._lookup_one("exomes", chrom, pos, ref, alt)
        gn_af = self._lookup_one("genomes", chrom, pos, ref, alt)

        out["exomes_af"] = ex_af
        out["genomes_af"] = gn_af

        vals = [x for x in (ex_af, gn_af) if isinstance(x, (int, float))]
        if vals:
            out["af"] = max(vals)
            out["source"] = "gnomAD-local"

        return out


# ══════════════════════════════════════════════════════════════
# ClinGen Dosage Sensitivity
# ══════════════════════════════════════════════════════════════

class ClinGenAnnotator:
    """ClinGen 유전자 용량 민감도 데이터 조회"""

    def __init__(self, clingen_tsv_path: str):
        self._data: Dict[str, Dict[str, Any]] = {}
        if clingen_tsv_path and os.path.exists(clingen_tsv_path):
            self._load(clingen_tsv_path)
        else:
            logger.warning(f"ClinGen TSV not found: {clingen_tsv_path}")

    def _load(self, path: str):
        def norm(x: str) -> str:
            return (x or "").strip()

        def to_int(x):
            try:
                return int(str(x).strip())
            except Exception:
                return None

        with open(path, "r", encoding="utf-8", errors="replace") as f:
            header_cols = None
            for line in f:
                if not line.strip():
                    continue
                if line.startswith("#Gene Symbol\t"):
                    header_cols = line.lstrip("#").rstrip("\n").split("\t")
                    break

            if not header_cols:
                logger.warning("ClinGen header not found")
                return

            reader = csv.DictReader(f, delimiter="\t", fieldnames=header_cols)
            for row in reader:
                gene = norm(row.get("Gene Symbol") or row.get("geneSymbol") or row.get("Gene") or "").upper()
                if not gene:
                    continue

                hi_score = to_int(row.get("Haploinsufficiency Score") or row.get("HI Score") or row.get("haploScore"))
                ts_score = to_int(row.get("Triplosensitivity Score") or row.get("TS Score") or row.get("triploScore"))

                candidate = {
                    "hgnc_id": norm(row.get("HGNC ID") or row.get("HGNCID") or ""),
                    "hi_score": hi_score,
                    "ts_score": ts_score,
                    "hi_desc": norm(row.get("Haploinsufficiency Description") or row.get("haploDescription") or ""),
                    "ts_desc": norm(row.get("Triplosensitivity Description") or row.get("triploDescription") or ""),
                    "last_eval": norm(row.get("Date Last Evaluated") or row.get("dateLastEvaluated") or ""),
                    "hi_disease_id": norm(row.get("Haploinsufficiency Disease ID") or ""),
                    "ts_disease_id": norm(row.get("Triplosensitivity Disease ID") or ""),
                    "url": norm(row.get("url") or row.get("URL") or ""),
                }

                existing = self._data.get(gene)
                if existing is None:
                    self._data[gene] = candidate
                else:
                    old_rank = (existing.get("hi_score") is not None, existing.get("ts_score") is not None)
                    new_rank = (candidate.get("hi_score") is not None, candidate.get("ts_score") is not None)
                    if new_rank > old_rank:
                        self._data[gene] = candidate

        logger.info(f"ClinGen loaded: {len(self._data)} genes")

    def lookup(self, gene: str) -> Dict[str, Any]:
        return self._data.get((gene or "").strip().upper(), {})


# ══════════════════════════════════════════════════════════════
# MANE RefSeq Map
# ══════════════════════════════════════════════════════════════

class MANEAnnotator:
    """MANE RefSeq 표준 전사체 매핑"""

    def __init__(self, mane_gff_path: str):
        self._data: Dict[str, Dict[str, str]] = {}
        if mane_gff_path and os.path.exists(mane_gff_path):
            self._load(mane_gff_path)
        else:
            logger.warning(f"MANE GFF not found: {mane_gff_path}")

    def _load(self, path: str):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "transcript":
                    continue

                attrs: Dict[str, str] = {}
                for kv in parts[8].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        attrs[k.strip()] = v.strip()

                gene = attrs.get("gene_name")
                enst = attrs.get("transcript_id")
                nm = ""
                dbx = attrs.get("Dbxref", "")
                if dbx:
                    for item in dbx.split(","):
                        item = item.strip()
                        if item.startswith("RefSeq:NM_"):
                            nm = item.replace("RefSeq:", "")
                            break

                if not gene:
                    continue
                entry = self._data.setdefault(gene, {})
                if enst and "enst" not in entry:
                    entry["enst"] = enst
                if nm and "nm" not in entry:
                    entry["nm"] = nm

        logger.info(f"MANE loaded: {len(self._data)} genes")

    def lookup(self, gene: str) -> Dict[str, str]:
        return self._data.get(gene, {})


# ══════════════════════════════════════════════════════════════
# Gene Intervals (BED-based positional gene lookup)
# ══════════════════════════════════════════════════════════════

class GeneIntervalAnnotator:
    """BED 기반 위치 → 유전자 매핑"""

    def __init__(self, gene_bed_path: str):
        self._trees: Dict[str, Any] = {}
        if gene_bed_path and os.path.exists(gene_bed_path):
            self._load(gene_bed_path)
        else:
            logger.warning(f"Gene BED not found: {gene_bed_path}")

    def _load(self, path: str):
        try:
            from intervaltree import IntervalTree
        except ImportError:
            logger.warning("intervaltree not installed, positional gene lookup disabled")
            return

        trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                chrom, start, end, gene = parts[0], int(parts[1]), int(parts[2]), parts[3].strip()
                trees[chrom].addi(start, end + 1, gene)

        self._trees = dict(trees)
        total = sum(len(t) for t in self._trees.values())
        logger.info(f"Gene intervals loaded: {total} intervals")

    def lookup(self, chrom: str, pos: int) -> Optional[str]:
        tree = self._trees.get(chrom)
        if not tree:
            return None
        hits = tree[pos]
        if hits:
            return next(iter(hits)).data
        return None


# ══════════════════════════════════════════════════════════════
# dbSNP URL Helper
# ══════════════════════════════════════════════════════════════

def make_dbsnp_url(rsid: str) -> str:
    """dbSNP rsID에서 NCBI URL을 생성합니다."""
    if rsid is None:
        return ""
    s = str(rsid).strip()
    if not s or s in (".", "0"):
        return ""
    m = re.search(r"\brs(\d+)\b", s, flags=re.IGNORECASE)
    if m:
        return f"https://www.ncbi.nlm.nih.gov/snp/?term={m.group(1)}"
    m2 = re.search(r"\b(\d+)\b", s)
    if m2:
        return f"https://www.ncbi.nlm.nih.gov/snp/?term={m2.group(1)}"
    return ""


# ══════════════════════════════════════════════════════════════
# Unified Annotator (Facade)
# ══════════════════════════════════════════════════════════════

class VariantAnnotator:
    """
    모든 annotation 소스를 통합하는 Facade 클래스.

    사용법:
        annotator = VariantAnnotator(config)
        result = annotator.annotate(chrom, pos, ref, alt, gene, ...)
    """

    def __init__(
        self,
        clinvar_vcf: str = "",
        gnomad_dir: str = "",
        gnomad_genomes_glob: str = "gnomad.genomes.v*.sites*.bgz",
        gnomad_exomes_glob: str = "gnomad.exomes.v*.sites*.bgz",
        clingen_tsv: str = "",
        mane_gff: str = "",
        gene_bed: str = "",
    ):
        self.clinvar = ClinVarAnnotator(clinvar_vcf)
        self.gnomad = GnomADAnnotator(gnomad_dir, gnomad_genomes_glob, gnomad_exomes_glob)
        self.clingen = ClinGenAnnotator(clingen_tsv)
        self.mane = MANEAnnotator(mane_gff)
        self.gene_intervals = GeneIntervalAnnotator(gene_bed)

        logger.info("VariantAnnotator initialized")

    def annotate(
        self,
        chrom: str, pos: int, ref: str, alt: str,
        gene: str = "", transcript: str = "",
        hgvsc: str = "", hgvsp: str = "", effect: str = "",
        sample_metrics: Dict[str, Any] = None,
    ) -> Dict[str, Any]:
        """
        단일 변이에 대해 모든 annotation을 수행합니다.

        Returns:
            phenotype_portal과 동일한 필드 구조의 딕셔너리
        """
        # MANE
        mane = self.mane.lookup(gene) if gene else {}

        # ClinGen
        clingen = self.clingen.lookup(gene) if gene else {}
        hgnc_id = clingen.get("hgnc_id", "")
        clingen_url = ""
        if hgnc_id:
            clingen_url = f"https://search.clinicalgenome.org/kb/genes/{hgnc_id.replace(':', '%3A')}"
        elif clingen.get("url"):
            clingen_url = clingen["url"]

        # gnomAD
        gnomad = self.gnomad.lookup(chrom, pos, ref, alt)
        gnomad_af = gnomad.get("af")

        # ClinVar
        clinvar = self.clinvar.lookup(chrom, pos, ref, alt)
        cv_primary = (clinvar or {}).get("clnsig_primary", "")

        # dbSNP rsID: ClinVar rsID 우선, 없으면 빈 문자열
        dbsnp_rsid = (clinvar or {}).get("rsid", "")
        dbsnp_url = make_dbsnp_url(dbsnp_rsid)

        # 샘플 메트릭스 기본값
        sm = sample_metrics or {}

        return {
            # 위치
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,

            # 유전자/전사체
            "gene": gene,
            "transcript": transcript,
            "canonical_enst": mane.get("enst", ""),
            "clinical_nm": mane.get("nm", ""),
            "hgvsc": hgvsc,
            "hgvsp": hgvsp,
            "effect": effect,

            # 샘플 메트릭스
            "dp": sm.get("dp"),
            "ref_depth": sm.get("ref_depth"),
            "alt_depth": sm.get("alt_depth"),
            "vaf": sm.get("vaf"),
            "gt": sm.get("gt", ""),
            "zygosity": sm.get("zygosity", ""),

            # gnomAD
            "gnomad_af": gnomad_af,
            "gnomad_exomes_af": gnomad.get("exomes_af"),
            "gnomad_genomes_af": gnomad.get("genomes_af"),
            "gnomad_source": gnomad.get("source", ""),

            # ClinVar
            "clinvar_sig": (clinvar or {}).get("clnsig", ""),
            "clinvar_sig_primary": cv_primary,
            "clinvar_conflicting": bool((clinvar or {}).get("conflicting", False)),
            "clinvar_conflict_detail": (clinvar or {}).get("conflict_summary", ""),
            "clinvar_revstat": (clinvar or {}).get("revstat", ""),
            "clinvar_stars": int((clinvar or {}).get("stars", 0) or 0),
            "clinvar_dn": (clinvar or {}).get("clndn", ""),
            "clinvar_variation_id": (clinvar or {}).get("variation_id", ""),

            # dbSNP
            "dbsnp_rsid": dbsnp_rsid,
            "dbsnp_url": dbsnp_url,

            # ClinGen
            "clingen_hi_score": clingen.get("hi_score"),
            "clingen_ts_score": clingen.get("ts_score"),
            "clingen_hi_desc": clingen.get("hi_desc", ""),
            "clingen_ts_desc": clingen.get("ts_desc", ""),
            "clingen_last_eval": clingen.get("last_eval", ""),
            "clingen_hi_disease_id": clingen.get("hi_disease_id", ""),
            "clingen_ts_disease_id": clingen.get("ts_disease_id", ""),
            "clingen_hgnc_id": hgnc_id,
            "clingen_url": clingen_url,
        }
