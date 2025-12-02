## MycoSNP Usage <a id="top"></a>

This guide explains how to prepare inputs and run both MycoSNP workflows: rapid taxonomic classification / subtyping (Pre-MycoSNP) and full reference-based variant analysis (NFCORE_MYCOSNP). It reflects repository state as of 2025-11-14.

### Table of Contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Updating & Reproducibility](#updating--reproducibility)
4. [Selecting a Workflow (`--mode`)](#selecting-a-workflow--mode)
5. [Pipeline Parameters](#parameters)
6. [Samplesheet Format](#samplesheet-format)
   - [Multiple Runs / Lanes](#multiple-runs--lanes)
   - [Integrated SRA Accessions](#integrated-sra)
   - [Integrating Existing VCFs](#integrated-vcfs)
   - [Automated Samplesheet Creation](#auto-samplesheet)
7. [Pre-MycoSNP-specific Inputs](#premycosnp-inputs)
8. [Reference Options (NFCORE_MYCOSNP)](#reference-options)
9. [Running Pre-MycoSNP](#run-premycosnp)
10. [Running NFCORE_MYCOSNP](#run-main)
11. [AMD-P QC Extensions](#amdp)
12. [Common Nextflow Arguments](#common-nextflow-arguments)
13. [Resource / Configuration Customization](#resource--configuration-customization)
14. [Background Execution](#background-execution)
15. [Memory Limits for Nextflow JVM](#memory-limits-for-nextflow-jvm)
16. [Deprecated Parameters](#deprecated)

---

### Requirements <a id="requirements"></a>

| Component         | Required Version / Notes                                                                             |
| ----------------- | ---------------------------------------------------------------------------------------------------- |
| Nextflow          | >= 24.10.4 (per `manifest.nextflowVersion`)                                                          |
| Java              | 17+                                                                                                  |
| Bash              | 3.2+                                                                                                 |
| Container runtime | Docker / Apptainer(Singularity) / Podman / etc. Conda or Mamba also supported but less reproducible. |
| RAM / CPU         | Adjust via profiles; high-depth assemblies and joint genotyping benefit from ≥8 CPUs / ≥32 GB RAM.   |

> [!TIP]
> If Apptainer/Singularity mount errors occur ("No such file or directory"), test an alternative runtime (Docker) or verify host kernel / FUSE settings.

### Installation <a id="installation"></a>

Install Nextflow first (see official docs). Then either clone or run directly:

```bash
git clone https://github.com/CDCgov/mycosnp-nf
```

Run from remote without cloning (cached under `~/.nextflow/assets/`):

```bash
nextflow run CDCgov/mycosnp-nf -profile singularity,test
```

### Updating & Reproducibility <a id="updating"></a>

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull CDCgov/mycosnp-nf
```

## Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [CDCgov/mycosnp-nf releases page](https://github.com/CDCgov/mycosnp-nf/releases) and find the latest version number (eg. `v1.6.4`). Then specify this when running the pipeline with `-r` (one hyphen) - e.g. `-r v1.6.4`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

### Selecting a Workflow (`--mode`) <a id="workflow-selection"></a>

The pipeline contains two DSL2 workflows. Switching is controlled by the `--mode` parameter (not `--workflow`, which is deprecated). Internally this drives inclusion / exclusion of reference-dependent variant calling steps.

- `PRE_MYCOSNP`: No reference required; produces taxonomic & subtype summaries plus assemblies and QC.
- `NFCORE_MYCOSNP`: Requires a reference or indices; produces full variant, consensus, phylogeny, optional annotation.
  If `--mode` is omitted it defaults to `NFCORE_MYCOSNP`.

Use `--mode PRE_MYCOSNP` for rapid taxon + subtype; omit or set `--mode NFCORE_MYCOSNP` for full variant analysis.

### Pipeline Parameters <a id="parameters"></a>

Quick help:

```bash
nextflow run CDCgov/mycosnp-nf --help
nextflow run CDCgov/mycosnp-nf --help --show_hidden
```

The following sections detail all user-facing parameters grouped by function.

#### Core Input / Output Parameters

| Param                           | Type       | Default                  | Notes                                                                                |
| ------------------------------- | ---------- | ------------------------ | ------------------------------------------------------------------------------------ |
| `--input`                       | file (csv) | (required)               | Samplesheet with FASTQ, SRA accessions, and/or VCF entries (multi-modal ingestion).  |
| `--outdir`                      | directory  | `./results`              | Root results directory. Absolute path recommended on cloud systems.                  |
| `--publish_dir_mode`            | string     | `copy`                   | Nextflow publish mode (`copy`, `link`, etc.). Hidden in help unless `--show_hidden`. |
| `--tracedir`                    | directory  | `<outdir>/pipeline_info` | Execution reports, timeline, trace, DAG. Auto-generated.                             |
| `--multiqc_title`               | string     | (null)                   | Title injected into MultiQC report.                                                  |
| `--multiqc_config`              | file       | (null)                   | Custom MultiQC YAML config. Hidden by default.                                       |
| `--multiqc_logo`                | file       | (null)                   | Logo referenced in MultiQC config. Hidden by default.                                |
| `--multiqc_methods_description` | file       | (null)                   | Methods HTML/YAML included in MultiQC.                                               |

#### Pre-MycoSNP Assembly & Classification Parameters

| Param              | Default | Description                                                                                                               |
| ------------------ | ------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------ |
| `--assembler`      | `skesa` | Shovill assembler (`skesa`, `spades`, `megahit`, `velvet`).                                                               |
| `--shovill_depth`  | 70      | Downsample target depth passed to Shovill.                                                                                |
| `--genome_size`    | ''      | Approx genome size; empty = auto / skip.                                                                                  |
| `--min_contig_cov` | 10      | Minimum contig coverage retained.                                                                                         |
| `--min_contig_len` | 300     | Minimum contig length retained.                                                                                           |
| `--gambit_db`      | null    | Path to GAMBIT metadata db (`.gdb`/`.db`). Required for taxonomic classification.                                         |
| `--gambit_h5_dir`  | null    | Directory of GAMBIT signature chunks (`\*.gs                                                                              | \*.h5[.gz]`). Required with `--gambit_db`. |
| `--subtype_db`     | null    | Directory containing `sourmash_taxa.csv` and signature files for subtype prediction. Optional but enables subtype output. |

#### Reference Genome & Masking Parameters

Provide either: (A) `--fasta` raw reference (masking + indexing performed) OR (B) `--ref_dir` pre-built directory OR (C) all of `--ref_masked_fasta --ref_fai --ref_bwa --ref_dict`.

| Param                | Default                      | Description                                                                                                                      |
| -------------------- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `--fasta`            | null                         | Raw reference FASTA; triggers masking (if `--mask true`) and index builds.                                                       |
| `--mask`             | true                         | Perform repeat masking (NUCMER + BEDTOOLS). Set false to skip.                                                                   |
| `--ref_dir`          | null                         | Use previously generated reference artifact directory (contains masked, fai, dict, bwa folder). Overrides individual ref params. |
| `--ref_masked_fasta` | null                         | Masked reference FASTA (used when bypassing build). Must accompany fai, bwa, dict.                                               |
| `--ref_fai`          | null                         | Samtools FASTA index.                                                                                                            |
| `--ref_bwa`          | null                         | BWA index directory (named `bwa`).                                                                                               |
| `--ref_dict`         | null                         | Picard sequence dictionary.                                                                                                      |
| `--genome`           | null                         | nf-core iGenomes key (alternative to manual reference specification). Ignored if `--fasta` or ref indices provided.              |
| `--igenomes_base`    | `s3://ngi-igenomes/igenomes` | Base path for iGenomes.                                                                                                          |
| `--igenomes_ignore`  | false                        | Skip loading iGenomes configuration.                                                                                             |

#### Variant Calling & Filtering Parameters

| Param                      | Default                                                                                           | Description                                                                                                                 |
| -------------------------- | ------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `--sample_ploidy`          | 1                                                                                                 | Ploidy for GATK HaplotypeCaller.                                                                                            |
| `--coverage`               | 0                                                                                                 | Target coverage for downsampling (0 = no downsampling).                                                                     |
| `--gvcfs_filter`           | `QD < 2.0 \|\| FS > 60.0 \|\| MQ < 40.0 \|\| DP < 10`                                             | GATK VariantFiltration expression.                                                                                          |
| `--gatkgenotypes_filter`   | `--min_GQ "50" --keep_GQ_0_refs --min_percent_alt_in_AD "0.8" --min_total_DP "10" --keep_all_ref` | Parameters for `filterGatkGenotypes.py`.                                                                                    |
| `--max_amb_samples`        | 10000000                                                                                          | Max ambiguous samples allowed.                                                                                              |
| `--max_perc_amb_samples`   | 10                                                                                                | Max % ambiguous samples allowed.                                                                                            |
| `--min_depth`              | 10                                                                                                | Minimum depth for consensus; bases below become `N` (affects QC report column header and phylogeny genome fraction metric). |
| `--skip_combined_analysis` | false                                                                                             | If true: stop after per-sample gVCF generation (no combined joint genotyping / phylogeny).                                  |

#### Phylogeny & Distance Matrix Parameters

| Param              | Default | Description                                                                            |
| ------------------ | ------- | -------------------------------------------------------------------------------------- |
| `--rapidnj`        | true    | Generate RapidNJ tree (neighbour-joining).                                             |
| `--fasttree`       | true    | Generate FastTree approximate maximum likelihood tree.                                 |
| `--iqtree`         | false   | Optional IQ-TREE build (slower, model selection).                                      |
| `--raxmlng`        | false   | Optional RAxML-NG build (slower, bootstrap).                                           |
| `--skip_phylogeny` | false   | Skip all tree-building (still produces SNP distance if not skipped combined analysis). |

#### Annotation Parameters (snpEff / snpeffr)

| Param            | Default                                                              | Description                                                              |
| ---------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------ |
| `--snpeff`       | false                                                                | Run snpEff annotation + snpeffr report generation.                       |
| `--snpeffconfig` | null                                                                 | Path to snpEff config directory.                                         |
| `--snpeffcache`  | null                                                                 | Path to snpEff cached database (required if `--snpeff true`).            |
| `--species`      | candida_auris_gca_016772135.1                                        | Species/id used for snpEff db selection & custom FKS1 hotspot reporting. |
| `--genes`        | CAB11_002014                                                         | Gene identifiers passed to snpeffr.                                      |
| `--positions`    | fks1_hs1=221637:221663,fks1_hs2=223782:223805,fks1_hs3=221805:221807 | Hotspot coordinates for FKS1 variant reporting.                          |
| `--exclude`      | synonymous_variant                                                   | Variant consequence types to ignore in FKS1 report.                      |

#### Sample Inclusion / Skipping Parameters

| Param                 | Default | Description                                                                    |
| --------------------- | ------- | ------------------------------------------------------------------------------ |
| `--skip_samples`      | ""      | Comma-separated sample IDs excluded from joint combine / downstream analyses.  |
| `--skip_samples_file` | null    | File with newline-separated sample IDs to skip (merged with `--skip_samples`). |

#### AMD-P (Public Health) Extension Parameters

Enabled by default (`--amdp true`). Adds QC parsing and threshold evaluation.

| Param               | Default                                                        | Description                                                                      |
| ------------------- | -------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| `--amdp`            | true                                                           | Activate AMD-P specific outputs (QC parsing, thresholds, metadata hooks).        |
| `--qc_thresholds`   | `GCrangePct:42-47.5,AvgQscore:28,RefLenCov:20,MeanCovDepth:20` | Key:value or range list evaluated by `qc_parser.py` to derive pass/fail summary. |
| `--metadata_csv`    | ""                                                             | Optional metadata file (future ingestion / reporting).                           |
| `--geolocation_csv` | null                                                           | Optional geolocation data for downstream integration.                            |
| `--test_samples`    | null                                                           | Restrict to listed samples for testing scenarios.                                |
| `--percent_n`       | 8                                                              | Threshold for percentage of Ns in consensus (affects QC evaluation).             |

#### Email & Reporting Parameters

| Param                      | Default | Description                                                       |
| -------------------------- | ------- | ----------------------------------------------------------------- |
| `--email`                  | null    | Address for completion summary email.                             |
| `--email_on_fail`          | null    | Email only on failed run.                                         |
| `--plaintext_email`        | false   | Send plain text instead of HTML email.                            |
| `--max_multiqc_email_size` | 25.MB   | Attach MultiQC if below this size. Hidden unless `--show_hidden`. |

#### Logging & Help Parameters

| Param               | Default | Description                                    |
| ------------------- | ------- | ---------------------------------------------- |
| `--help`            | false   | Print standard parameter help.                 |
| `--help_full`       | false   | Extended nf-schema help (shows hidden params). |
| `--show_hidden`     | false   | Display hidden parameters in help output.      |
| `--validate_params` | true    | Enforce schema validation at launch.           |
| `--monochrome_logs` | false   | Disable ANSI color in logs.                    |
| `--version`         | false   | Print pipeline version and exit.               |

#### Institutional Config / Advanced Parameters

Hidden parameters typically set via profiles:

| Param                            | Default                                                  | Purpose                                   |
| -------------------------------- | -------------------------------------------------------- | ----------------------------------------- |
| `--custom_config_version`        | master                                                   | nf-core institutional config commit.      |
| `--custom_config_base`           | https://raw.githubusercontent.com/nf-core/configs/master | Base URL for institutional config fetch.  |
| `--config_profile_name`          | null                                                     | Profile name annotation.                  |
| `--config_profile_description`   | null                                                     | Profile description annotation.           |
| `--config_profile_contact`       | null                                                     | Profile contact annotation.               |
| `--config_profile_url`           | null                                                     | Profile URL annotation.                   |
| `--pipelines_testdata_base_path` | https://raw.githubusercontent.com/nf-core/test-datasets/ | Base for test datasets.                   |
| `--trace_report_suffix`          | timestamp                                                | Added to report filenames for uniqueness. |

#### Execution Profiles

Profiles defined in `nextflow.config` tailor execution environment:

| Profile                                                               | Key Overrides                                                               |
| --------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `premycosnp`                                                          | `mode=PRE_MYCOSNP`                                                          |
| `singleSample`                                                        | `coverage=70`, `skip_phylogeny=true`, `snpeff=true`, `igenomes_ignore=true` |
| `vcfs`                                                                | `amdp=true`, `igenomes_ignore=true`                                         |
| `debug`                                                               | Enables hash dump, disables cleanup, turns on process name validation.      |
| `docker` / `singularity` / `conda` / `mamba` / `podman` / `apptainer` | Container / package manager activation.                                     |
| `wave`                                                                | Enables Wave frozen strategy for environments.                              |
| `gitpod`                                                              | Local resource limits for Gitpod workspace.                                 |
| `gpu`                                                                 | GPU runtime flags.                                                          |

#### Parameter Interactions & Validation

| Scenario                                                                            | Behavior                                                                                             |
| ----------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| Provide `--ref_dir`                                                                 | Ignores `--fasta` and all `--ref_*` individual index params.                                         |
| Provide any of `--ref_masked_fasta --ref_fai --ref_bwa --ref_dict` (no `--ref_dir`) | Must supply all four; otherwise run exits with error.                                                |
| `--skip_combined_analysis true`                                                     | Stops after per-sample gVCFs; no joint genotyping, filtering, phylogeny, or annotation.              |
| `--skip_phylogeny true`                                                             | Suppresses tree building (RapidNJ/FastTree/IQ-TREE/RAxML-NG); distance matrix may still be produced. |
| `--snpeff true` without `--snpeffcache`                                             | Will fail; cache is required.                                                                        |
| `--amdp false`                                                                      | QC thresholds not evaluated; `combined_QC_results.csv` not produced.                                 |
| `--coverage 0`                                                                      | Disables downsampling (SEQTK step skipped).                                                          |
| `--mask false`                                                                      | Skips repeat masking; reference used as-is.                                                          |

#### Troubleshooting Tips

| Issue                      | Check                                                                        |
| -------------------------- | ---------------------------------------------------------------------------- |
| Missing phylogeny outputs  | Was `--skip_phylogeny` set or trees disabled (`--rapidnj/--fasttree false`)? |
| Empty subtype CSVs         | Species rank not determined or taxon absent from `sourmash_taxa.csv`.        |
| No annotation outputs      | `--snpeff` not set or missing `--snpeffcache`.                               |
| Reference build skipped    | Presence of `--ref_dir` overrides build from `--fasta`.                      |
| QC combined results absent | `--amdp` disabled or thresholds file misformatted.                           |

### Samplesheet Format <a id="samplesheet-format"></a>

Single integrated samplesheet handles FASTQ, SRA accessions, and existing VCF inclusion. Schema: `assets/schema_input.json`.

Minimum requirement: `sample` plus one of `fastq_1`, `sra`, or `vcf`.

Columns:
| Column | Required? | Description |
|--------|-----------|-------------|
| `sample` | yes | Unique sample ID (no spaces; converted internally if present). |
| `fastq_1` | conditional | R1 gzipped FASTQ (`.fastq.gz` / `.fq.gz`). |
| `fastq_2` | conditional | R2 gzipped FASTQ; presence triggers paired-end. |
| `sra` | conditional | SRA accession ID; pipeline downloads and ingests reads. |
| `vcf` | conditional | Existing per-sample gVCF/VCF (`.vcf.gz`) produced with the SAME reference used in current run. Index (`.tbi`) must be adjacent. |

Example mixed samplesheet (see updated comprehensive example at `assets/samplesheet.csv`):

```csv
sample,fastq_1,fastq_2,sra,vcf
CAURIS_01,/data/reads/CAURIS_01_R1.fastq.gz,/data/reads/CAURIS_01_R2.fastq.gz,,
CAURIS_02,, ,SRR1234567,
LEGACY_03,,,,/data/prev_runs/LEGACY_03.g.vcf.gz
```

#### Multiple Runs / Lanes <a id="multiple-lanes"></a>

Provide additional lane FASTQs as extra columns (`fastq_3`, `fastq_4`, ...) with consistent sample ID; pipeline merges prior to QC.

#### Integrated SRA Accessions <a id="integrated-sra"></a>

Place the accession in the `sra` column. No separate `--add_sra_file` is used (deprecated). Sample ID will label downloaded FASTQs.

#### Integrating Existing VCFs <a id="integrated-vcfs"></a>

Specify each VCF path in the `vcf` column. The run will include those in combined genotyping / downstream phylogeny unless `--skip_combined_analysis true`. Former `--add_vcf_file` is deprecated.

### Pre-MycoSNP-specific Inputs <a id="premycosnp-inputs"></a>

All Pre-MycoSNP inputs are read-level (FASTQ or SRA) plus optional subtyping database. No reference is needed because GAMBIT classification operates on assembly sketches. Ensure both `--gambit_db` and `--gambit_h5_dir` are supplied; failure to provide both results in an early exit.

- GAMBIT database: `--gambit_db` path to metadata `.gdb/.db` file.
- GAMBIT signatures: `--gambit_h5_dir` directory containing split signature chunks (`*.gs|*.h5[.gz]`).
- Subtyping database: `--subtype_db` directory with `sourmash_taxa.csv` + signature files; `candida_auris_clades.sig` supports clade differentiation.

### Reference Options (NFCORE_MYCOSNP) <a id="reference-options"></a>

Supply ONE of the following strategies:

1. `--fasta reference.fa` (raw FASTA; masking + index build performed if `--mask true`).
2. `--ref_dir path/to/previous/reference/` (contains subfolders: `masked/`, `dict/`, `fai/`, `bwa/`).
3. All four explicit files: `--ref_masked_fasta`, `--ref_fai`, `--ref_bwa`, `--ref_dict`.

If any explicit ref index is given, all must be present or the run exits.

Optional: `--genome` iGenomes key (ignored if other reference inputs supplied).

### Running Pre-MycoSNP <a id="run-premycosnp"></a>

```bash
nextflow run CDCgov/mycosnp-nf --mode PRE_MYCOSNP \
  --input samplesheet.csv \
  --gambit_db gambit.db --gambit_h5_dir signatures/ \
  --subtype_db sourmash_db/ \
  -profile singularity \
  --outdir results_premycosnp
```

Key Outputs: assemblies (`samples/<id>/assembly/`), taxonomy CSVs, subtype CSVs, per-sample line summaries, combined summary, MultiQC.

### Running NFCORE_MYCOSNP <a id="run-main"></a>

Minimum viable run requires: samplesheet with at least one read source (FASTQ or SRA) and `--fasta` OR pre-built indices. For large cohorts consider providing pre-built indices via `--ref_dir` to save setup time.

Basic full run example (variant analysis + phylogeny):

```bash
nextflow run CDCgov/mycosnp-nf \
  --input samplesheet.csv \
  --fasta reference.fa \
  --snpeff true --snpeffcache snpeff_cache/ \
  --iqtree true --raxmlng true \
  -profile singularity \
  --outdir results_full
```

Test profiles:

```bash
nextflow run CDCgov/mycosnp-nf -profile test,singularity
nextflow run CDCgov/mycosnp-nf -profile test_full,singularity
```

Work directory layout:

```text
work/          # Intermediate process files (hash-keyed)
results/       # Published outputs (controlled by --outdir)
.nextflow.log  # Nextflow execution log
```

### AMD-P QC Extensions <a id="amdp"></a>

Enabled by default (`--amdp true`): parses QC into `combined_QC_results.csv` using thresholds in `--qc_thresholds`. Disable with `--amdp false` for leaner output.

### Common Nextflow Arguments <a id="nextflow-args"></a>

`-profile`, `-resume`, and `-c` behave per standard nf-core pipelines. Multiple profiles may be comma-separated (`-profile test,docker`). If no profile is specified, local execution expects tools on `PATH`.

### Resource / Configuration Customization <a id="resource-customization"></a>

Common labels (see `conf/base.config`): `process_low`, `process_medium`, `process_high`. Override by process name or label group.
Example scaling all high processes:

```nextflow
process {
  withLabel: process_high { cpus = 16; memory = 64.GB }
}
```

Dynamic retry logic multiplies resources for certain exit codes (e.g. 137 OOM) up to 3 attempts before halting.
Increase memory / CPUs via a custom config using `withName:` selectors or labels. Example to boost FAQCS memory:

```nextflow
process {
  withName: FAQCS { memory = 50.GB }
}
```

If a process repeatedly exits with code 137 (OOM), raise `memory`; for 143 (timeout) adjust `time`.

### Background Execution <a id="background"></a>

Use `-bg`, `screen`, or `tmux`. Retain `work/` if you intend to resume or reuse intermediate artifacts for differential analysis.

### Memory Limits for Nextflow JVM <a id="jvm-memory"></a>

Set JVM caps to avoid runaway allocation:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

### Deprecated Parameters <a id="deprecated"></a>

| Deprecated       | Replacement                 | Rationale                                     |
| ---------------- | --------------------------- | --------------------------------------------- |
| `--workflow`     | `--mode`                    | Unified naming with config param key.         |
| `--add_sra_file` | `sra` column in samplesheet | Simplifies ingestion, single source of truth. |
| `--add_vcf_file` | `vcf` column in samplesheet | Consolidates external VCF inclusion.          |

---

## General nf-core documentation

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](/conf/base.config#L17) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the mycosnp-nf pipeline is failing after multiple re-submissions of the `BWA_PRE_PROCESS:FAQCS` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `MYCOSNP:BWA_PRE_PROCESS:FAQCS (SAMPLE_1` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'MYCOSNP:BWA_PRE_PROCESS:FAQCS (SAMPLE_1'

Caused by:
    Process `MYCOSNP:BWA_PRE_PROCESS:FAQCS (SAMPLE_1)` terminated with an error exit status (137)

Command executed:
FaQCs \
            -d . \
            -u SAMPLE_1.fastq.gz \
            --prefix SAMPLE_1_clean \
            -t 8 \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    FaQCs -d . -u SAMPLE_1.fastq.gz -t 8 <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [CDCgov/mycosnp-nf Github repo](https://github.com/CDCgov/mycosnp-nf/search?q=process+FAQCS). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/modules/faqcs/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_medium`](/modules/nf-core/modules/faqcs/main.nf#L3). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organize workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_medium` label are set in the pipeline's [`base.config`](/conf/base.config#L32-L36) which in this case is defined as 16GB. Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `FAQCS` process failure by creating a custom config file that sets at least 16GB of memory, in this case increased to 50GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: FAQCS {
        memory = 50.GB
    }
}
```

> **NB:** We specify just the process name i.e. `FAQCS` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `MYCOSNP:BWA_PRE_PROCESS:FAQCS`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)

2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)

3. Create the custom config accordingly:
   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [nf-core Slack](https://nf-co.re/join/slack) in the `#configs` channel.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
export NXF_OPTS='-Xms1g -Xmx4g'
```
