# CDCgov/mycosnp-nf: Output

## Introduction

This guide describes the outputs produced by the two selectable modes of the pipeline:

1. PRE_MYCOSNP (rapid taxonomic classification + subtype + light QC)
2. NFCORE_MYCOSNP (reference-based variant calling, QC, and phylogeny; optional annotation)

The pipeline publishes files into a structured directory tree rooted at `<OUTDIR>` (parameter `--outdir`). Below is a high-level map (subdirectories omitted for brevity). Actual directories are created only when the corresponding process runs (e.g. phylogeny trees skipped if `--skip_phylogeny`).

### NFCORE_MYCOSNP Output Structure (Reference-based variant calling mode)

```text
<OUTDIR>/
├── aggregate_outputs/                # Cross-sample aggregated reports
│   ├── combined_QC_results.csv       # QC threshold evaluation results (if --amdp true)
│   └── combined_cauris_refB11205_fks1.csv  # snpeffr FKS1 hotspot report (if --snpeff true)
├── combined/                         # Combined (multi-sample) variant calling and phylogeny outputs
│   ├── filteredgvcfs/                # Filtered combined genotype VCF files with indices
│   ├── finalfiltered/                # Final filtered VCF for downstream analysis
│   ├── genotypegvcfs/                # Combined genotyped gVCF files with indices
│   ├── gvcf/                         # Combined gVCF files with indices
│   ├── microreact/                   # Microreact metadata CSV for visualization (if generated)
│   ├── phylogeny/                    # Tree files organized by method
│   │   ├── fasttree/                 # FastTree phylogeny outputs (if --fasttree true)
│   │   ├── iqtree/                   # IQ-TREE phylogeny outputs (if --iqtree true)
│   │   ├── quicksnp/                 # QuickSNP tree outputs (if generated)
│   │   ├── raxmlng/                  # RAxML-NG phylogeny outputs (if --raxmlng true)
│   │   └── rapidnj/                  # RapidNJ phylogeny outputs (if --rapidnj true)
│   ├── qc_report/                    # Combined QC report from QC_REPORTSHEET
│   ├── samtools_coverage/            # Per-sample coverage summary files
│   ├── samtools_depth/               # Per-sample depth files
│   ├── samtools_flagstat/            # Per-sample flagstat files
│   ├── samtools_idxstats/            # Per-sample index statistics
│   ├── samtools_stats/               # Per-sample comprehensive alignment statistics
│   ├── selectedsnps/                 # Selected SNP VCF files with indices
│   ├── selectedsnpsfiltered/         # Filtered selected SNPs with indices
│   ├── snpdists/                     # SNP distance matrix (combined.tsv)
│   ├── snpeff/                       # snpEff annotation outputs (if --snpeff true)
│   ├── splitvcf/                     # Per-sample split VCF files and sample list
│   ├── vcf-qc-report/                # VCF-level QC summary and passed samples list
│   └── vcf-to-fasta/                 # SNP alignment FASTA (vcf-to-fasta.fasta)
├── multiqc/                          # Unified MultiQC report
│   ├── multiqc_data/                 # MultiQC data files and logs
│   ├── multiqc_plots/                # Plot exports (pdf, png, svg subdirectories)
│   └── multiqc_report.html           # Interactive HTML report
├── pipeline_info/                    # Nextflow execution reports and metadata
│   ├── execution_report_*.html       # Execution report with resource usage
│   ├── execution_timeline_*.html     # Timeline visualization
│   ├── execution_trace_*.txt         # Detailed process execution trace
│   ├── params_*.json                 # Run parameters
│   ├── pipeline_dag_*.html           # Pipeline DAG visualization
│   └── software_versions.yml         # Software versions used
├── reference/                        # Masked reference + indices (if --save_reference true)
│   ├── bwa/                          # BWA index files
│   ├── dict/                         # Picard sequence dictionary
│   ├── fai/                          # Samtools FASTA index
│   └── masked/                       # Masked reference FASTA and masking coordinates
└── samples/                          # Per-sample outputs
    └── <sample_id>/
        ├── faqcs/                    # Trimmed reads and trimming statistics
        ├── fastqc_post/              # FastQC reports on trimmed reads
        ├── fastqc_raw/               # FastQC reports on raw reads
        ├── finalbam/                 # Final coordinate-sorted BAM with index
        ├── qc_report/                # Per-sample QC summary
        └── variant_calling/          # Per-sample variant calling outputs
            └── haplotypecaller/      # GATK HaplotypeCaller gVCF files with indices
```

### PRE_MYCOSNP Output Structure (Rapid taxonomic classification mode)

```text
<OUTDIR>/
├── aggregate_outputs/                # Cross-sample aggregated reports
│   └── pre-mycosnp_summary/
│       └── pre-mycosnp-summary.csv   # Combined summary of all samples
├── multiqc/                          # Unified MultiQC report
│   ├── multiqc_data/                 # MultiQC data files and logs
│   ├── multiqc_plots/                # Plot exports (pdf, png, svg subdirectories)
│   └── multiqc_report.html           # Interactive HTML report
├── pipeline_info/                    # Nextflow execution reports and metadata
│   ├── execution_report_*.html       # Execution report with resource usage
│   ├── execution_timeline_*.html     # Timeline visualization
│   ├── execution_trace_*.txt         # Detailed process execution trace
│   ├── params_*.json                 # Run parameters
│   ├── pipeline_dag_*.html           # Pipeline DAG visualization
│   └── software_versions.yml         # Software versions used
└── samples/                          # Per-sample outputs
    └── <sample_id>/
        ├── assembly/                 # De novo assembly FASTA
        ├── faqcs/                    # Trimmed reads and trimming statistics
        ├── fastqc_raw/               # FastQC reports on raw reads
        ├── pre-mycosnp_summary/      # Per-sample line summary CSV
        ├── subtype/                  # Subtype prediction CSV (if applicable)
        └── taxonomy/                 # GAMBIT taxonomic classification CSV
```

Key conditional outputs:

| Condition / Param          | Affected Outputs                                                                                     |
| -------------------------- | ---------------------------------------------------------------------------------------------------- |
| `--mode PRE_MYCOSNP`       | Pre-MycoSNP sample assembly, taxonomy, subtype, pre-mycosnp summaries (no reference-based alignment) |
| `--skip_combined_analysis` | Skips combining gVCFs and downstream variant consolidation (no combined/\* variant directories)      |
| `--skip_phylogeny`         | Skips tree-building subdirectories under `combined/phylogeny/`                                       |
| `--snpeff`                 | Generates per-sample snpEff outputs + snpeffr aggregate CSV                                          |
| `--amdp` (default false)   | Runs QC_PARSER producing `qc/combined_QC_results.csv` (enable with `--amdp true`)                    |
| `--save_reference false`   | Suppresses publishing of `reference/*` indices                                                       |
| `--save_alignment false`   | Suppresses publication of trimmed reads, final BAMs, alignment QC outputs                            |
| `--save_debug true`        | Publishes additional intermediate/debug BAMs and stats (picard*\* directories, seqkit*\* etc.)       |

## Navigation

- [Pre-MycoSNP Mode](#pre-mycosnp-mode)
  - [Sample-level QC](#sample-level-qc-pre)
  - [Assembly](#assembly)
  - [Taxonomic classification (GAMBIT)](#taxonomic-classification-gambit)
  - [Subtyping (sourmash)](#subtyping-sourmash)
  - [Pre-MycoSNP combined summary](#pre-mycosnp-combined-summary)
- [NFCORE_MYCOSNP](#nfcore_mycosnp)
  - [Reference preparation](#reference-preparation)
  - [Alignment + preprocessing](#alignment--preprocessing)
  - [Combined QC report](#combined-qc-report)
  - [Variant calling and consolidation](#variant-calling-and-consolidation)
  - [VCF QC + SNP FASTA + Distances](#vcf-qc--snp-fasta--distances)
  - [Phylogeny](#phylogeny)
  - [Annotation (optional snpEff)](#annotation-optional-snpeff)
  - [MultiQC](#multiqc)
  - [Pipeline information](#pipeline-information)

---

## Pre-MycoSNP Mode

> [!NOTE]
> Pre-MycoSNP taxonomic classification and _Candida auris_ clade typing have been validated for common fungal pathogens. See the validation report on the Wiki.

### Sample-level QC (Pre)

Tools: FastQC (raw), SeqKit (pair filtering), FaQCs (trimming + statistics).

Per-sample directories under `samples/<sample_id>/`:

| Directory                   | Contents                                                                                                                                            |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `fastqc_raw/`               | Raw FastQC reports (`*.html`, `*.zip`)                                                                                                              |
| `faqcs/`                    | Trimmed reads (`*_1.trimmed.fastq.gz`, `*_2.trimmed.fastq.gz`), stats (`*.stats.txt`), per-base quality files (`*.for_qual_histogram.txt` if debug) |
| `seqkit_pair/` (debug only) | Unpaired filtering intermediate fastqs                                                                                                              |

### Assembly

Shovill performs de novo assembly (default assembler SKESA; depth capped by `--shovill_depth` = 70; correction skipped via `--nocorr`).

| Path                            | File                                         |
| ------------------------------- | -------------------------------------------- |
| `samples/<sample_id>/assembly/` | `<sample_id>.fa` (renamed from `contigs.fa`) |

### Taxonomic classification (GAMBIT)

GAMBIT sketch querying against provided fungal database chunks (`--gambit_db` + `--gambit_h5_dir`). Distances reported as Jaccard distance (1 - Jaccard index). Species prediction attempted; falls back to genus; else reports closest match only.

| Path                            | File                     |
| ------------------------------- | ------------------------ |
| `samples/<sample_id>/taxonomy/` | `<sample_id>_gambit.csv` |

### Subtyping (sourmash)

Performed only when predicted taxon present and matching entries exist in `assets/sourmash_db/sourmash_taxa.csv` with corresponding signature(s) in `assets/sourmash_db/signatures/`.

For _Candida auris_, clade (I–VI) estimated using `candida_auris_clades.sig`; ANI estimated via `sourmash compare --containment --ani`.

| Path                           | File                                                   |
| ------------------------------ | ------------------------------------------------------ |
| `samples/<sample_id>/subtype/` | `<sample_id>_subtype.csv` (may be empty if no subtype) |

### Pre-MycoSNP combined summary

All per-sample line summaries (`*_linesummary.csv`) from `samples/<sample_id>/pre-mycosnp_summary/` are combined.

| Path                                     | File                      |
| ---------------------------------------- | ------------------------- |
| `aggregate_outputs/pre-mycosnp_summary/` | `pre-mycosnp-summary.csv` |

Header columns (current implementation – all prefixed with `PM_` except `Sample`):

`Sample,PM_Predicted_Rank,PM_Predicted_Taxon,PM_Subtype_Closest_Match,PM_Subtype_ANI,PM_Closest_GAMBIT_Entry_Description,PM_Closest_GAMBIT_Entry_Distance,PM_Trimmed_Reads,PM_Avg_Read_Quality,PM_Sample_Assembly_Length,PM_Sample_Assembly_GC,PM_Reference_Genome_Length,PM_Avg_Depth_Coverage,PM_Reference_GC`

Changes vs older documentation: all data columns now have a `PM_` prefix; paths moved from `combined/pre-mycosnp_summary` to `aggregate_outputs/pre-mycosnp_summary`.

> [!NOTE]
> For _Candida auris_, if `PM_Subtype_ANI < 99.7`, `PM_Subtype_Closest_Match` is replaced with the caution message about the clade separation threshold (99.7) per validation study.

Column descriptions:

| Column                              | Description                                                   |
| ----------------------------------- | ------------------------------------------------------------- |
| Sample                              | Sample identifier                                             |
| PM_Predicted_Rank                   | GAMBIT predicted rank: `species`, `genus`, or `no prediction` |
| PM_Predicted_Taxon                  | GAMBIT predicted taxon (blank if none)                        |
| PM_Subtype_Closest_Match            | sourmash closest subtype match or caution message             |
| PM_Subtype_ANI                      | Estimated ANI (%) to subtype signature                        |
| PM_Closest_GAMBIT_Entry_Description | Description of closest GAMBIT database entry                  |
| PM_Closest_GAMBIT_Entry_Distance    | Jaccard distance (1 - Jaccard index)                          |
| PM_Trimmed_Reads                    | Reads after trimming (FaQCs)                                  |
| PM_Avg_Read_Quality                 | Average Q score post trimming                                 |
| PM_Sample_Assembly_Length           | Total assembly length (bp)                                    |
| PM_Sample_Assembly_GC               | Assembly GC (%)                                               |
| PM_Reference_Genome_Length          | Length of closest reference genome (if species predicted)     |
| PM_Avg_Depth_Coverage               | Total trimmed bases / reference genome length (approx depth)  |
| PM_Reference_GC                     | GC (%) of closest reference                                   |

---

## NFCORE_MYCOSNP

### Reference preparation

If `--fasta` (or `--ref_dir` / explicit indices) provided, repeats are masked (`nucmer` + BED masking) and indices generated:

| Directory           | Files                                                                                 |
| ------------------- | ------------------------------------------------------------------------------------- |
| `reference/masked/` | `reference.fa` (masked FASTA) + optional debug `*.coords` / `*.bed` if `--save_debug` |
| `reference/bwa/`    | `bwa` index directory containing `reference.{amb,ann,bwt,pac,sa}`                     |
| `reference/fai/`    | `reference.fa.fai`                                                                    |
| `reference/dict/`   | `reference.dict`                                                                      |

### Alignment + preprocessing

Per-sample alignment and preprocessing steps (lane merge, pair filtering, optional downsampling if `--coverage > 0`, trimming, alignment, duplicate marking, cleaning, mate fixing, read group assignment, indexing). Key publication directories (some only with `--save_debug true`):

| Directory                                            | Purpose                                                  |
| ---------------------------------------------------- | -------------------------------------------------------- |
| `samples/<sample_id>/fastqc_raw/`                    | Raw FastQC                                               |
| `samples/<sample_id>/faqcs/`                         | Trimmed reads + stats                                    |
| `samples/<sample_id>/bwamem/` (debug)                | Initial aligned BAM (unsorted or intermediate)           |
| `samples/<sample_id>/picard_markduplicates/` (debug) | Duplicate-marked BAM                                     |
| `samples/<sample_id>/picard_cleansam/` (debug)       | Cleaned BAM                                              |
| `samples/<sample_id>/picard_fixmate/` (debug)        | Mate-fixed BAM                                           |
| `samples/<sample_id>/finalbam/`                      | Final coordinate-sorted, RG added BAM + `.bai`           |
| `samples/<sample_id>/fastqc_post/`                   | FastQC on filtered/trimmed reads                         |
| `combined/samtools_stats/`                           | `*.stats` summaries (comprehensive alignment statistics) |
| `combined/samtools_flagstat/`                        | `*.flagstat` files (read flag statistics)                |
| `combined/samtools_idxstats/`                        | `*.idxstats` files (index statistics)                    |
| `combined/samtools_coverage/`                        | `*.txt` coverage summaries (mean depth, base quality)    |
| `combined/samtools_depth/`                           | `*.tsv` per-position depth files                         |

### Combined QC report

Generated by `QC_REPORTSHEET` from FaQCs trimming statistics and Samtools alignment metrics (`samtools coverage`, `samtools depth`, and `samtools stats`). Path:

| Path                  | File            |
| --------------------- | --------------- |
| `combined/qc_report/` | `qc_report.txt` |

If `--amdp true` is enabled, the parsed threshold evaluation appears at:

| Path  | File                                                                      |
| ----- | ------------------------------------------------------------------------- |
| `qc/` | `combined_QC_results.csv` (pass/fail per thresholds in `--qc_thresholds`) |

QC Report columns (dynamic depth column label reflects `--min_depth`):

`Sample Name, Reads Before Trimming, GC Before Trimming, Average Q Score Before Trimming, Reference Length Coverage Before Trimming, Reads After Trimming, Paired Reads After Trimming, Unpaired Reads After Trimming, GC After Trimming, Average Q Score After Trimming, Reference Length Coverage After Trimming, Mean Coverage Depth, Reads Mapped, Genome Fraction at <min_depth>X`

**Metric Sources:**

- Trimming metrics (reads, GC, Q scores): `FaQCs`
- Mean Coverage Depth: `samtools coverage` (weighted mean across all references)
- Reads Mapped: `samtools stats` (count and percentage)
- Genome Fraction at depth threshold: `samtools depth` (positions ≥ min_depth / total reference length)

### Variant calling and consolidation

If not skipped (`--skip_combined_analysis false`):

1. Per-sample gVCFs: `samples/<sample_id>/variant_calling/haplotypecaller/` (`*.g.vcf.gz` + index)
2. Combined gVCF: `combined/gvcf/` (`combined.combined.g.vcf.gz`)
3. Genotyping: `combined/genotypegvcfs/` (`combined_genotype.vcf.gz`)
4. Filtering: `combined/filteredgvcfs/` (`combined_genotype_filtered.vcf.gz`) using expression from `--gvcfs_filter`
5. SNP selection: `combined/selectedsnps/` (`combined_genotype_filtered_snps.vcf.gz`)
6. Genotype-level filtering: `combined/selectedsnpsfiltered/` (`combined_genotype_filtered_snps_filtered.vcf.gz`) applying thresholds from `--gatkgenotypes_filter`
7. Final compressed view: `combined/finalfiltered/` (`finalfiltered.vcf.gz` + index)
8. Individual sample splits (SNPs): `<OUTDIR>/<sample_id>/splitvcf/`

Additional outputs:

| Directory                     | Files                                            |
| ----------------------------- | ------------------------------------------------ |
| `combined/vcf-qc-report/`     | `*.report.txt`, `*.samples.txt` (VCF QC summary) |
| `combined/consensus/` (debug) | `*.fasta.gz` consensus sequences (optional)      |

### VCF QC + SNP FASTA + Distances

| Directory                | Purpose                                                                        |
| ------------------------ | ------------------------------------------------------------------------------ |
| `combined/vcf-to-fasta/` | Alignment FASTA (`vcf-to-fasta.fasta`) sites after SNP selection and filtering |
| `combined/snpdists/`     | SNP pairwise distance matrix (`*.tsv`)                                         |

### Phylogeny

Built only if `--skip_phylogeny false`. Depending on enabled methods (defaults: RapidNJ + FastTree):

| Directory                      | Tree File (standardized renamed)                                  |
| ------------------------------ | ----------------------------------------------------------------- |
| `combined/phylogeny/rapidnj/`  | `rapidnj_phylogeny.nh`                                            |
| `combined/phylogeny/fasttree/` | `fasttree_phylogeny.nh`                                           |
| `combined/phylogeny/iqtree/`   | `iqtree_phylogeny.nh` (if `--iqtree true`)                        |
| `combined/phylogeny/raxmlng/`  | `raxmlng_bestTree.nh`, `raxmlng_support.nh` (if `--raxmlng true`) |
| `combined/phylogeny/quicksnp/` | `quicksnp_phylogeny.nwk` (if quicksnp enabled)                    |

### Annotation (optional snpEff)

Enabled with `--snpeff true` (requires cache/config parameters). Per-sample:

| Directory                     | Files                                                                           |
| ----------------------------- | ------------------------------------------------------------------------------- |
| `samples/<sample_id>/snpeff/` | snpEff outputs (`*.csv`, `*.txt`, `*.html`, `<sample_id>.snpeff.vcf.gz`, index) |
| `aggregate_outputs/`          | `snpeffr` FKS1 hotspot report (`*.csv`)                                         |

> [!NOTE]
> The `snpeffr` report may show `undetermined` mutation entries where variants disappear between per-sample gVCF and combined genotype stages.

### MultiQC

Path: `multiqc/`

Files:

| File                                                                | Description                                 |
| ------------------------------------------------------------------- | ------------------------------------------- |
| `multiqc_report.html`                                               | Standalone aggregated QC report             |
| `multiqc_data/`                                                     | Parsed metric JSON / TSV from modules       |
| `multiqc_plots/`                                                    | Exported static plot images                 |
| `workflow_summary_mqc.yaml`                                         | Parameter summary injected for traceability |
| `premycosnp_software_versions.yml`, `mycosnp_software_versions.yml` | Included in report data set if present      |

FastQC plots here reflect raw (untrimmed) reads; trimming improvements appear in FaQCs statistics.

### Pipeline information

Path: `pipeline_info/`

| File                                  | Purpose                                           |
| ------------------------------------- | ------------------------------------------------- |
| `execution_report_<timestamp>.html`   | Run summary report                                |
| `execution_timeline_<timestamp>.html` | Process timeline visualization                    |
| `execution_trace_<timestamp>.txt`     | Tabular trace (resources, statuses)               |
| `pipeline_dag_<timestamp>.html`       | DAG graph of process dependencies                 |
| `premycosnp_software_versions.yml`    | Collated tool versions (Pre mode)                 |
| `mycosnp_software_versions.yml`       | Collated tool versions (NFCORE mode)              |
| `samplesheet.valid.csv`               | Sanitized / validated samplesheet used internally |
| `params.json`                         | Parameters snapshot                               |

---

## GAMBIT & sourmash Databases (Assets)

Database manifests reside under `assets/` (e.g. `assets/gambit_db/` and `assets/sourmash_db/`). Ensure you provide `--gambit_db` and `--gambit_h5_dir` pointing to the GAMBIT CSV + H5 fragment directory, and `--subtype_db` pointing to sourmash subtype signatures. The workflow reconstructs GAMBIT H5 file (`gambit_signatures.gs`) transiently and removes it after querying.

---

## Troubleshooting & Tips

| Symptom                                     | Likely Cause                                    | Action                                                         |
| ------------------------------------------- | ----------------------------------------------- | -------------------------------------------------------------- |
| Missing `combined/gvcf` directory           | `--skip_combined_analysis` set                  | Omit the flag to perform joint genotyping                      |
| Empty subtype CSV                           | Taxon not in sourmash DB or rank not species    | Confirm entry in `sourmash_taxa.csv` and signature presence    |
| No phylogeny trees                          | `--skip_phylogeny true` or no SNPs              | Disable skip flag; verify SNP FASTA not empty                  |
| Low `PM_Subtype_ANI` (<99.7) for _C. auris_ | Genuine ambiguity or sequencing/assembly issues | Re-sequence or review read quality; caution interpreting clade |
| Consensus FASTA missing                     | Debug outputs disabled                          | Set `--save_debug true` if you need consensus sequences        |

---

## Citation

Please cite the MycoSNP pipeline (CDCgov/mycosnp-nf) and underlying tools (GAMBIT, sourmash, Shovill, BWA, GATK, snpEff, MultiQC) as listed in `CITATIONS.md`.

---

## Change Summary (Relative to Previous Documentation)

1. Updated directory structure (added `aggregate_outputs/`, granular `combined/*` folders, removed monolithic `stats/`).
2. Corrected Pre-MycoSNP summary path and column names (`PM_` prefixes).
3. Adjusted QC report path from `stats/qc_report/` to `combined/qc_report/` plus optional parsed thresholds in `qc/`.
4. Added VCF QC report, consensus FASTA, snpdists, and per-method phylogeny directories.
5. Clarified conditional publication controlled by parameters (`--snpeff`, `--skip_phylogeny`, etc.).
6. Differentiated raw vs trimmed FastQC context and integrated workflow parameter summary file.
7. Explicit mapping for variant consolidation step outputs.

---

## Glossary

| Term             | Definition                                        |
| ---------------- | ------------------------------------------------- |
| gVCF             | Genomic VCF with per-site reference confidence    |
| ANI              | Average Nucleotide Identity percentage            |
| Jaccard distance | 1 - Jaccard index (k-mer sketch distance metric)  |
| Consensus FASTA  | Per-sample sequence incorporating called variants |

---

## Feedback

Open issues or feature requests via the repository issue templates. Include relevant paths (e.g. `combined/qc_report/qc_report.txt`) and parameter JSON (`pipeline_info/params.json`) for reproducibility.
