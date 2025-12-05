## MycoSNP Documentation

MycoSNP is a portable Nextflow workflow for whole-genome sequencing analysis of fungal pathogens (with a current focus on _Candida auris_). It supports two complementary modes:

1. **Pre-MycoSNP (`--mode PRE_MYCOSNP`)** – Fast taxonomic classification (GAMBIT) and clade / subtype prediction (sourmash) plus minimal assembly/QC summaries.
2. **NFCORE_MYCOSNP (`--mode NFCORE_MYCOSNP`)** – Full reference-based variant calling (GATK), filtering, consensus construction, phylogeny, optional annotation (snpEff + snpeffr), and AMD-P public health QC extensions.

This README gives a high-level orientation. See linked pages for full details.

### Table of Contents

1. [Key Features](#key-features)
2. [Documentation Pages](#documentation-pages)
3. [Workflows Overview](#workflows-overview)
4. [Prerequisites](#prerequisites)
5. [Installation](#installation)
6. [Quick Start Examples](#quick-start-examples)
7. [AMD-P Extensions](#amd-p-extensions)
8. [Typical Outputs](#typical-outputs)
9. [Parameter Reference](#parameter-reference)
10. [Troubleshooting](#troubleshooting)
11. [Citation & Attribution](#citation--attribution)
12. [Support & Contributing](#support--contributing)

---

### Key Features

- Dual-mode analysis: rapid classification vs. full variant workflow.
- Automated repeat masking & indexing for new references (NUCMER, BEDTOOLS, BWA, Picard, Samtools).
- gVCF generation + joint variant combining and SNP-only filtering.
- Multiple phylogeny methods (RapidNJ, FastTree, optional IQ-TREE / RAxML-NG).
- Subtyping via sourmash signatures (e.g., _C. auris_ clade prediction with ANI threshold logic).
- AMD-P QC threshold parsing and combined pass/fail reports.
- Optional snpEff annotation and FKS1 hotspot reporting via snpeffr.
- Consolidated MultiQC report including workflow summary metadata.
- Reproducible container / conda environment support across multiple execution profiles.

### Documentation Pages

| Page                        | Description                                                                                     |
| --------------------------- | ----------------------------------------------------------------------------------------------- |
| [`usage.md`](usage.md)      | Running the pipeline, samplesheet format, execution profiles, and complete parameter reference. |
| [`output.md`](output.md)    | Detailed explanation of directory structure and generated files.                                |
| `dev_notes.md` (if present) | Developer guidelines and maintenance notes.                                                     |

### Workflows Overview

#### PRE_MYCOSNP Mode

**Purpose:** Rapid species identification and subtype prediction

**Core Steps:**

1. FastQC (raw reads)
2. SeqKit pairing
3. FaQCs (quality control and trimming)
4. Shovill assembly
5. GAMBIT taxonomic classification
6. Subtype prediction (sourmash)
7. Summary statistics
8. MultiQC aggregation

**Principal Outputs:**

- Per-sample taxon summary CSVs
- Combined Pre-MycoSNP summary report
- Taxonomy and subtype classification tables
- Draft assemblies

#### NFCORE_MYCOSNP

**Purpose:** Comprehensive variant calling and phylogenetic analysis

**Core Steps:**

1. Reference preparation (repeat masking and indexing)
2. FastQC (raw reads)
3. Read quality control and trimming
4. Read alignment and BAM quality assessment
5. Per-sample gVCF generation (GATK HaplotypeCaller)
6. Joint variant calling and filtering
7. Consensus sequence and SNP FASTA generation
8. Phylogenetic tree construction (RapidNJ/FastTree)
9. Optional variant annotation (snpEff)
10. QC reporting and MultiQC aggregation

**Principal Outputs:**

- Per-sample gVCF files
- Combined filtered SNP VCF
- Consensus FASTA sequences
- SNP distance matrix
- Phylogenetic trees (Newick format)
- QC reports (pass/fail thresholds)
- Variant annotation outputs (if enabled)
- MultiQC report with software versions

### Prerequisites

- **Nextflow** ≥ version 24.10.4.
- One supported execution profile: Docker, Singularity/Apptainer, Podman, Conda/Mamba, or Wave.
- Adequate RAM / CPU for assembly and variant calling (adjust process resources via profiles).
- Reference FASTA (unless using iGenomes key or supplying pre-built indices via `--ref_dir`).
- GAMBIT DB + signature chunks for Pre-MycoSNP classification (`--gambit_db`, `--gambit_h5_dir`).
- sourmash subtype database (optional for clade/subtype prediction).
- snpEff cache (required if `--snpeff true`).

### Installation

Clone the repository or rely on hosted revision:

```bash
git clone https://github.com/CDCgov/mycosnp-nf
```

Pull containers automatically (Docker / Singularity) when running Nextflow.

### Quick Start Examples

Pre-MycoSNP only:

```bash
nextflow run CDCgov/mycosnp-nf -profile singularity \
  --input samplesheet.csv --mode PRE_MYCOSNP \
  --gambit_db gambit.db --gambit_h5_dir signatures/ \
  --subtype_db sourmash_db/ --outdir results_premycosnp
```

Full workflow with phylogeny and annotation:

```bash
nextflow run CDCgov/mycosnp-nf -profile singularity \
  --input samplesheet.csv --fasta reference.fa \
  --snpeff true --snpeffcache snpeff_cache/ \
  --iqtree true --raxmlng true --outdir results_full
```

Generate per-sample gVCFs only (no joint analysis):

```bash
nextflow run CDCgov/mycosnp-nf -profile singularity \
  --input samplesheet.csv --fasta reference.fa \
  --skip_combined_analysis true --outdir results_gvcf_only
```

### AMD-P Extensions

Disabled by default. Enable with `--amdp true` to add:

- QC threshold evaluation (`combined_QC_results.csv`).
- Metadata and geolocation hooks (future integration).
- Percent N and coverage thresholds into pass/fail logic.

### Typical Outputs (High-Level)

See [`output.md`](output.md) for full tree. Key locations:

- `samples/<sample_id>/faqcs/` – Trimmed read QC.
- `samples/<sample_id>/finalbam/` – Final BAM + index (NFCORE_MYCOSNP).
- `combined/gvcf/`, `combined/finalfiltered/`, `combined/phylogeny/*/` – Variant and phylogeny artifacts.
- `aggregate_outputs/pre-mycosnp_summary/` – Combined Pre-MycoSNP summary.
- `aggregate_outputs/` – Annotation aggregate outputs (e.g. FKS1 report) when snpEff enabled.
- `pipeline_info/` – Execution reports, trace, timeline, DAG, software versions YAML.
- `multiqc/` – Aggregated QC and workflow summary.

### Parameter Reference

See [`usage.md`](usage.md) for detailed parameter tables, groupings, defaults, and interactions.

### Troubleshooting

| Symptom                  | Suggestion                                                                                           |
| ------------------------ | ---------------------------------------------------------------------------------------------------- |
| Missing phylogeny trees  | Verify `--skip_phylogeny` not set and specific tree flags (`--rapidnj`, `--fasttree`, etc.) enabled. |
| Empty subtype CSV        | Species rank not resolved or taxon absent from `sourmash_taxa.csv`.                                  |
| Annotation files absent  | Ensure `--snpeff true` and `--snpeffcache` path supplied.                                            |
| Reference not rebuilt    | `--ref_dir` supplied; remove to trigger masking/index.                                               |
| Combined QC file missing | `--amdp` disabled or malformed `--qc_thresholds`.                                                    |

### Citation & Attribution

Please cite MycoSNP (CDCgov/mycosnp-nf) and underlying tools (GAMBIT, sourmash, GATK, RapidNJ, FastTree, IQ-TREE, RAxML-NG, snpEff) in publications. Include the pipeline version (see MultiQC or `mycosnp_software_versions.yml`). A formal citation entry will be added when a DOI is minted.

### Support & Contributing

- Open issues for bugs or feature requests on the repository issue tracker.
- Follow coding style and module patterns; see nf-core DSL2 guidelines.
- Run tests (`nf-test`) where provided before submitting PRs.
- Security or vulnerability concerns: use the dedicated PR template or disclosure channels.

---

For broader nf-core usage, installation, and configuration guidance visit the [nf-core website](https://nf-co.re).

Last updated: 2025-11-14
