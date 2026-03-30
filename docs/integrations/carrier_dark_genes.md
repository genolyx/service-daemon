# Carrier screening + unified dark-gene Nextflow

The **dark_gene_pipeline** repository is a **single** Nextflow workflow: FASTQ → align/dedup → panel SNV calling, CNV/SV, depth on `dark_genes_plus`, Paraphase/SMACA, ExpansionHunter, summary reports, etc. (see that repo’s `main.nf`).

The service-daemon **does not** vendor that pipeline. It can **drive** it directly when:

1. **`run_analysis.sh` is not used** — either remove/rename `src/run_analysis.sh` under the carrier layout, **or** set **`CARRIER_SCREENING_FORCE_DIRECT_NEXTFLOW=true`** so the daemon skips the wrapper even if that file exists. (Unset **`CARRIER_SCREENING_RUN_SCRIPT`** alone is not enough: the plugin auto-discovers `run_analysis.sh` under `CARRIER_SCREENING_PIPELINE_DIR` / the FASTQ layout parent.)
2. **`DARK_GENE_PIPELINE_DIR`** points at your checkout (directory containing `main.nf` and `nextflow.config`), **or** set **`CARRIER_SCREENING_MAIN_NF`** to that `main.nf`.
3. **Refs** are set in `.env` (`REF_FASTA`, `REF_FAI`, `REF_DICT`, optional `REF_BWA_INDICES`) so the daemon can pass `--ref_fasta`, `--ref_fai`, `--ref_dict`, etc.

## What the daemon passes

For this entry point the plugin uses **dark-gene parameter names**:

| Parameter | Source |
|-----------|--------|
| `fastq_dir` | Parent directory of `fastq_r1` (must contain paired FASTQs matching the pipeline’s glob) |
| `outdir` | Job `analysis_dir` |
| `ref_fasta`, `ref_fai`, `ref_dict`, `ref_bwa_indices` | Daemon settings |
| `backbone_bed` | Job / `CARRIER_DEFAULT_BACKBONE_BED` / `data/bed` resolution |
| `eh_catalog` | `<DARK_GENE_PIPELINE_DIR>/bed/variant_catalog_grch38.json` if present |
| `skip_cnv` | `DARK_GENE_SKIP_CNV` (default `true`; set `false` to run gCNV) |

Optional overrides from **`job.params`**: `pon_tar`, `gcnv_model`, `interval_list`, `hba_bed`, `cyp21a2_bed`, `eh_catalog`, `cleanup`.

Use **`CARRIER_SCREENING_NEXTFLOW_PROFILE`** (e.g. `docker`) if your pipeline requires `-profile`.

## Alternative: explicit main.nf

You can set **`CARRIER_SCREENING_MAIN_NF=/path/to/dark_gene_pipeline/main.nf`** instead of **`DARK_GENE_PIPELINE_DIR`**. The daemon still detects “dark gene” style by path substring `dark_gene` and uses the same CLI mapping.

## VCF for Review Case

The daemon’s **`check_completion`** searches for `*.vcf` / `*.vcf.gz` under the analysis tree. The dark pipeline publishes panel calls under **`outdir/variant/*_filtered.vcf.gz`**. Those are picked up for **annotation → `result.json`** like any other carrier run.

## SMAca C/T counts in the portal

The Review **SMAca CHECK** block shows **SNP C,T counts** when the detailed report (or merged sidecar) contains lines the portal parses, e.g. **`C_T=40,44`** or **`CT_counts=…`** (see `portal/index.html` / `dark_genes._smaca_extract_snp_ct_counts`).

**This repo does not run SMAca.** The **unified Nextflow** pipeline must **emit** those counts. Two supported options:

1. **Inline in `*_detailed_report.txt`** — In the SMAca section body, print:
   - `C_T=<C_reads>,<T_reads>` (or `CT_counts=…` / `SMN_C_T=…` per existing parsers).

2. **Sidecar file (recommended if SMAca already computes counts in code)** — Publish a small file next to the summary so the daemon can merge it into the detailed text before parsing:
   - **Path (first match wins):** `summary/smaca_snp_counts.txt`, `smaca/smaca_snp_counts.txt`, or any `**/smaca_snp_counts.txt` under the analysis outdir.
   - **Content (one line):** `C_T=40,44` or `40 44` (two integers).

   `service-daemon` reads this file during `collect_dark_genes_from_analysis_dir`, injects the line after `C_Ratio=` (or before `SMN1_CN=` if needed), and sets **`dark_genes.smaca_ct_sidecar`** to the relative path for audit.

After adding the sidecar or editing the detailed report, run **Reprocess only** (or a full re-run) so **`result.json`** is rebuilt.

## Review tab + PDF

- **`result.json`** includes **`dark_genes`**: populated when `summary/*_summary_report.txt` (or detailed) exists under **`analysis_dir`** after `process_results`.
- **Portal → Review → Dark genes** tab shows the summary/detailed text.
- **Generate Report** merges **`dark_genes`** from **`result.json`** into **`report.json`** as **`dark_genes.report_summary`** / **`report_detailed`** (truncated). Jinja templates **`carrier_EN.html`** / **`carrier_couples_EN.html`** add a **“Supplementary analysis (hard-to-sequence regions)”** page when `report_summary` is non-empty.

## When `run_analysis.sh` is used

If **`CARRIER_SCREENING_RUN_SCRIPT`** is set, the daemon **only** runs that wrapper and **does not** apply `DARK_GENE_PIPELINE_DIR` resolution. Use **either** the wrapper **or** direct Nextflow + `DARK_GENE_PIPELINE_DIR`, not both.
