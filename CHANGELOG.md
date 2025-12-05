# CDCgov/mycosnp-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.6.4 - [Unreleased]

This release represents the AMD-P (Advanced Molecular Detection Platform) enhancement of the CDCgov/mycosnp-nf pipeline with public health-specific enhancements and workflow improvements.

### `Added`

- **AMD-P QC Extensions**
  - `QC_PARSER` module for automated QC threshold evaluation
  - `qc_parser.py` script parses QC reports and generates pass/fail determinations based on configurable thresholds
  - `--qc_thresholds` parameter for customizable QC criteria (default: `GCrangePct:42-47.5,AvgQscore:28,RefLenCov:20,MeanCovDepth:20`)
  - `aggregate_outputs/combined_QC_results.csv` output with per-sample pass/fail status
  - AMD-P parameters: `--metadata_csv`, `--geolocation_csv`, `--test_samples`, `--percent_n`, `--amdp` (disabled by default, enable with `--amdp true`)

- **Microreact Integration**
  - `MICROREACTSHAPES` module for phylogenetic visualization data preparation
  - `update_microreact_shapes.py` script processes metadata and geolocation data for Microreact compatibility
  - Test sample filtering and geolocation enrichment support

- **Parameter Validation**
  - Migrated from nf-validation to nf-schema@2.4.2 plugin
  - Configurable parameter validation with ignore lists for non-standard parameters
  - Help system parameters: `--help`, `--help_full`, `--show_hidden`

- **Workflow Mode Selection**
  - `--mode` parameter for workflow selection
  - Supported modes: `PRE_MYCOSNP` and `NFCORE_MYCOSNP`

- **nf-test Coverage**
  - Unit tests added for all local modules
  - Pipeline-level tests added for `PRE_MYCOSNP` workflow
  - Profile-specific tests added for `singleSample` and `vcfs` execution profiles

### `Changed`

- **Parameter Validation Configuration**
  - Updated validation block to nf-schema-compliant syntax
  - Changed default `params.mode` to `NFCORE_MYCOSNP` (previously `params.workflow`)
  - Set `params.igenomes_ignore = true` by default

- **Samplesheet Format and Validation**
  - Unified samplesheet format consolidates FASTQ, SRA, and VCF inputs into single CSV
  - Required columns: `sample` plus at least one of `fastq_1`, `sra`, or `vcf`
  - Optional columns: `fastq_2` (paired-end)
  - Multi-lane support maintained (`fastq_3`, `fastq_4` are added as an additional row under column headers `fastq_1` and `fastq_2`, keeping the sample name the same)
  - SRA accessions integrated directly via `sra` column (eliminates separate `--add_sra_file` parameter)
  - Existing VCF files integrated via `vcf` column (eliminates separate `--add_vcf_file` parameter)
  - `INPUT_CHECK` subworkflow provides comprehensive validation with file existence checking
  - Automatic lane merging for samples with multiple FASTQ columns
  - Enhanced error reporting for missing or invalid input files

- **Module Migration to nf-core**
  - Replaced local modules with nf-core versions for improved maintainability:
    - `FAQCS` from nf-core/modules
    - `BWA/MEM` from nf-core/modules
    - `PICARD MARKDUPLICATES and ADDORREPLACEREADGROUPS` from nf-core/modules
    - `SEQTK/SAMPLE` from nf-core/modules (replaced local `seqtk_sample.nf`)
    - `GATK4/COMBINEGVCFS` from nf-core/modules (replaced local `gatk4_localcombinegvcfs.nf`)
  - Added `GATK4/INDEXFEATUREFILE` module from nf-core/modules

- **QC Metrics Transition**
  - Replaced Qualimap with Samtools tools for alignment quality metrics:
    - `SAMTOOLS/COVERAGE` for coverage statistics and mean depth calculations
    - `SAMTOOLS/DEPTH` for per-position depth analysis and genome fraction metrics
    - `SAMTOOLS/STATS` for comprehensive alignment statistics including mapped read counts
  - Updated `QC_REPORT` module to process Samtools outputs instead of Qualimap results
  - Maintained all QC report metrics while improving computational efficiency and consistency

- **Module Organization**
  - Reorganized local modules from flat structure (25 .nf files) to organized subdirectories (23 modules with main.nf/meta.yml structure)
  - Added `INPUT_PROC` local module for reference genome preprocessing
  - Updated module publication paths to `aggregate_outputs/` for cross-sample reports

- **Output Directory Structure**
  - QC parser outputs published to `aggregate_outputs/` directory
  - Consolidated cross-sample reports in single output location

### `Removed`

- **Local Modules Removed**
  - `seqtk_sample.nf`
  - `gatk4_localcombinegvcfs.nf`
  - `pre_mycosnp_comb_summary.nf`
  - `snpeff_ann.nf`, `snpeff_build.nf`, `snpeff_local.nf`
  - `faqcs.nf`
  - `bwa_mem.nf`
  - `picard_markduplicates.nf`
  - `picard_addorreplacereadgroups.nf`
  - `qualimap.nf` (replaced by Samtools tools for QC metrics)

- **Deprecated Perl and Python Scripts**
  - `check_samplesheet.py` (validation handled by nf-schema and INPUT_CHECK subworkflow)
  - `mycosnp_combine_lanes.pl` (lane merging implemented in INPUT_CHECK subworkflow)
  - `mycosnp_create_sample_sheet.pl`
  - `mycosnp_full_samplesheet.sh`
  - `fastq_dir_to_samplesheet.py`

- **Deprecated Parameters**
  - `--add_sra_file` parameter removed (SRA accessions now specified in `sra` column of samplesheet)
  - `--add_vcf_file` parameter removed (VCF files now specified in `vcf` column of samplesheet)

### `Deprecated`

- `--workflow` parameter deprecated in favor of `--mode` (backward compatibility maintained)
- nf-validation plugin replaced with nf-schema (no user action required)

### `Notes`

- AMD-P extensions (disabled by default, enable with `--amdp true`)
- Upstream features and fixes from v1.6.2 and earlier included
- Requires Nextflow >= 24.10.4 (updated from >= 21.10.3)
- Compatible with nf-schema@2.4.2 plugin
- Optimized for AWS Batch execution with AMDP profiles

### `Testing`

- Existing test profiles maintained and functional
- QC parser threshold parsing validated
- Documentation examples updated with current parameter names

---

## v1.6.3 - [2025-07-25]

### Fixed

- [Pre-MycoSNP workflow] GC content calculation for closest GAMBIT match (`Reference_GC` field in pre-mycosnp-summary.csv) excluded lowercase (soft-masked) bases, i.e. "c" and "g", leading to underestimates for `Reference_GC`. This fix includes lowercase "g" and "c" in the GC content calculation in `bin/pre-mycosnp-stats.sh`, for both the closest GAMBIT match and the sample assembly. However, the sample assembly shouldn't contain any soft-masked bases to begin with, so this shouldn't result in any differences in `Sample_Assembly_GC`.

---

## v1.6.2 - [2025-05-18]

### `Fixed`

- Docker permissions issue with GAMBIT container in some Linux environments (issue [#137](https://github.com/CDCgov/mycosnp-nf/issues/137))

### `Changed`

- Bumped GAMBIT fungal database from v0.2.0 to v1.0.0. See [GAMBIT'S documentation](https://theiagen.notion.site/GAMBIT-7c1376b861d0486abfbc316480046bdc#3f6610c81fbb4812b745234441514e12) for details about the databases.
- Bumped GAMBIT version from v1.0.0 to v1.1.0. Database is now included with the repo in `assets/gambit_db/`. Signatures database in `assets/gambit_db/signatures/` (split into chunks so each chunks is below GitHub's file size limit of 100 MB and warning limit of 50 MB). See GAMBIT v1.1.0 [release notes](https://github.com/jlumpe/gambit/releases/tag/v1.1.0) for more info (no functional changes).

## v1.6.1 - [2025-05-09]

### `Changed`

- Switched Shovill container image to StaPH-B's Docker image (`quay.io/staphb/shovill:1.1.0-2022Dec`). Addresses issues with pulling image from Docker Hub (#134).
- Reordered params in `nextflow_schema.json` to be more intuitive.

## v1.6.0 - [2025-04-21]

### `Added`

- Added the "Pre-MycoSNP" workflow for quick fungal taxonomic classification and _Candida auris_ clade typing (using de novo assemblies)
  - Base workflow development by Jared Johnson [@DOH-JDJ0303](https://github.com/DOH-JDJ0303) (#119)
  - Original sourmash subtyping subworkflow by Charlotte Royer [@royercj](https://github.com/royercj) (#114)
  - Co-development, validation, testing, and documentation by Zack Mudge [@zmudge3](https://github.com/zmudge3)
- Added the "Genome Fraction" metric to the main MycoSNP workflow's [QC report](/docs/output.md#qc-report-statsqc_reportqc_reporttxt), by CJ Jossart [@cjjossart](https://github.com/cjjossart) (#123, #128)
- Nextflow Tower / Seqera Cloud functionality, by Jared Johnson [@DOH-JDJ0303](https://github.com/DOH-JDJ0303) (#104)
- Added FKS1 coordinates for a third hotspot region (single amino acid) associated with echinocandin resistance in _Candida auris_ (for the snpeffr report - `results/snpeff/combined_cauris_refB11205_fks1.csv`)

### `Changed`

- Bumped GATK version to 4.5.0.0
- Changed SnpEff-related config to work with cloud, by Jared Johnson [@DOH-JDJ0303](https://github.com/DOH-JDJ0303) (#103)
- Updated samplesheet and reference genome links in [test.config](/conf/test.config)
- Changed filename of `results/snpeff/combined.csv` to `results/snpeff/combined_cauris_refB11205_fks1.csv`
- Bumped [snpeffr](https://github.com/CDCgov/snpeffr) version to v1.1.1. The "mutation" column in the snpeffr report (`results/snpeff/combined_cauris_refB11205_fks1.csv`) can now contain a value of "undetermined". This occurs when GATK identifies a variant in the individual VCF file for an isolate, but then the identified variant is lost in the subsequent step when GATK creates a merged VCF file that includes all isolates.
  - Snpeffr v1.1.1 updates by Mal Rajeev [@mrajeev08](https://github.com/mrajeev08)

### `Removed`

- Removed `--tmpdir` parameter and associated config in [nextflow.config](/nextflow.config)

### `Fixed`

- Fixed _Candida auris_ FKS1 hotspot 1 coordinates in [nextflow.config](/nextflow.config), to include the first nucleotide of F635, and remove two extra nucleotides at the end

### `Deprecated`

- Nextflow `-entry` parameter to run specific workflow is deprecated. Use MycoSNP's `--workflow` parameter instead.

## v1.5 Wingardium Leviosa - [2023-05-09]

- Re-released on 12/12/2024 following #126 (added empty tmp directory - no changes in functionality)
- Re-released on 09/03/2024 following #122 (disabled snpEff logging - no changes in functionality)

### `Fixed`

- Added new subworkflow for `snpEff` variant annotation
- Added `Quicksnp` module which creates a neighbor-joining tree with SNP distances on branches

### `Fixed`

- Corrected an issue with downsampling read count calculations

### `Deprecated`

- Removed `rate` parameter to use `coverage` instead

### `TODO`

- Update file naming for snpEff mutations file

## v1.4 Tremella Snidget - [2022-06-27]

### `Added`

- Added containers for processes which did not have any container specified.

### `Fixed`

- Changed default --min-depth to 10

## v1.3 Musky Albus - [2022-06-09]

### `Fixed`

- Changed downsample strategy in `modules/local/downsample_rate.nf` that was causing differences with results from geneflow version. Downsample rate now set at default 1 (closes #67)

## v1.2 Expecto Patronum - [2022-05-11]

### `Added`

- Changed minimum time requirements in `base.config` for `process_low`, `process_medium`, `process_high` to **72.h** and `process_long` to **120.h**.
- SNP distance matrix addition as output.
- Updated `qc_report.txt` to include coverage **mean depth** and **reads mapped**.
- Positions masked `(N)` based on **DP** & Added functionality to use `min_depth` (Default 50).
- Change `test` profile to include `min_depth = 2` so it will run to completion.

### `Fixed`

- Bug fix for downsample mismatch.
- Change configuration variable name from vcftools_filter to `gatkgenotypes_filter`.
- Changed samplesheet creation to accept multiple directories as arguments and to recursively search for sequences.
- Set full vcf consensus file to debug output
- Remove part nf-core branding

### `TODO`

- Update logo

## v1.1 Candid Aura - [2022-04-01]

### `Fixed`

- Fixed bug in `modules/local/lane_merge.nf` that was causing samplesheet CSV file to not recognize R2 (closes #39)
- Formatting of `docs/output.md`
- Changed output file `combined/vcf-to-fasta/combined_vcf-to-fasta.fasta` -> `combined/vcf-to-fasta/vcf-to-fasta.fasta`
- Output file `combined/vcf-to-fasta/vcf-to-fasta.fasta` will now replace stars `*` with dashes `-`
- Output file `combined/phylogeny/rapidnj/rapidnj_phylogeny.tre` -> `combined/phylogeny/rapidnj/rapidnj_phylogeny.nh`
- Output file `combined/phylogeny/iqtree/vcf-to-fasta.fasta.treefile` -> `combined/phylogeny/iqtree/iqtree_phylogeny.nh`
- Output file `combined/phylogeny/raxmlng/output.raxml.bestTree` -> `combined/phylogeny/raxmlng/raxmlng_bestTree.nh`
- Output file `combined/phylogeny/raxmlng/output.raxml.support` -> `combined/phylogeny/raxmlng/raxmlng_support.nh`

### `Deprecated`

- `/results/qc` output dir removed

### `TODO`

- Continue improving output docs

## v1.0 Espresso Myconaut - [2022-03-25]

Initial release of CDCgov/mycosnp-nf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added `test` and `test_full` profiles for easy testing
- Support for phylogenetic tree generation `--rapidnj` `--fasttree` `--iqtree` `--raxmlng`
- Skip reference generation `--ref_dir` | `assets/precomputed/reference` available for testing
- Add SRA ids to download from NCBI `--add_sra_file` | `assets/sra_small.csv` `assets/sra_large.csv` available for testing
- Add vcf files generated from previous runs `--add_vcf_file` | `assets/precomputed/vcf` available for testing
- Skip combined variant analysis (run reference prep and mapping) `skip_combined_analysis`
- Skip samples capability `--skip_samples`
- Skip samples file capability `--skip_samples_file`

### `Dependencies`

- Nextflow
- nf-core

### `Deprecated`

- GeneFlow support
- QC pdf report

### `TODO`

- Intermediate file cleanup and management
- Update logo and metro-style workflow
