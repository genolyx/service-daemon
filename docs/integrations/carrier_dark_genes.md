# Carrier screening + unified dark-gene Nextflow

The **dark_gene_pipeline** repository is a **single** Nextflow workflow: FASTQ → align/dedup → panel SNV calling, CNV/SV, depth on `dark_genes_plus`, Paraphase/SMACA, ExpansionHunter, summary reports, etc. (see that repo’s `main.nf`).

The service-daemon **does not** vendor that pipeline. It can **drive** it directly when:

1. **`CARRIER_SCREENING_RUN_SCRIPT`** is **unset** (no `run_analysis.sh` wrapper).
2. **`DARK_GENE_PIPELINE_DIR`** points at your checkout (directory containing `main.nf` and `nextflow.config`).
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

## Review tab + PDF

- **`result.json`** includes **`dark_genes`**: populated when `summary/*_summary_report.txt` (or detailed) exists under **`analysis_dir`** after `process_results`.
- **Portal → Review → Dark genes** tab shows the summary/detailed text.
- **Generate Report** merges **`dark_genes`** from **`result.json`** into **`report.json`** as **`dark_genes.report_summary`** / **`report_detailed`** (truncated). Jinja templates **`carrier_EN.html`** / **`carrier_couples_EN.html`** add a **“Supplementary analysis (hard-to-sequence regions)”** page when `report_summary` is non-empty.

## When `run_analysis.sh` is used

If **`CARRIER_SCREENING_RUN_SCRIPT`** is set, the daemon **only** runs that wrapper and **does not** apply `DARK_GENE_PIPELINE_DIR` resolution. Use **either** the wrapper **or** direct Nextflow + `DARK_GENE_PIPELINE_DIR`, not both.
