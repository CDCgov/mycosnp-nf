# CDCgov/mycosnp-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
