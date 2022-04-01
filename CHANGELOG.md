# nf-core/mycosnp: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1 Candid Aura - [04/01/2022]

### `Added`

### `Fixed`

*   Fixed bug in `modules/local/lane_merge.nf` that was causing samplesheet CSV file to not recognize R2 (closes #39)
*   Formatting of `docs/output.md`
*   Changed output file `combined/vcf-to-fasta/combined_vcf-to-fasta.fasta` -> `combined/vcf-to-fasta/vcf-to-fasta.fasta`
*   Output file `combined/vcf-to-fasta/vcf-to-fasta.fasta` will now replace stars `*` with dashes `-`
*   Output file `combined/phylogeny/rapidnj/rapidnj_phylogeny.tre` -> `combined/phylogeny/rapidnj/rapidnj_phylogeny.nh`
*   Output file `combined/phylogeny/iqtree/vcf-to-fasta.fasta.treefile` -> `combined/phylogeny/iqtree/iqtree_phylogeny.nh`
*   Output file `combined/phylogeny/raxmlng/output.raxml.bestTree` -> `combined/phylogeny/raxmlng/raxmlng_bestTree.nh`
*   Output file `combined/phylogeny/raxmlng/output.raxml.support` -> `combined/phylogeny/raxmlng/raxmlng_support.nh`

### `Dependencies`

### `Deprecated`

*   `/results/qc` output dir removed

### `TODO`

*	Continue improving output docs

---
## v1.0 Espresso Myconaut - [03/25/2022]

Initial release of CDCgov/mycosnp-nf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

*   Added `test` and `test_full` profiles for easy testing
*   Support for phylogenetic tree generation `--rapidnj` `--fasttree` `--iqtree` `--raxmlng`
*   Skip reference generation `--ref_dir` | `assets/precomputed/reference` available for testing
*   Add SRA ids to download from NCBI `--add_sra_file` | `assets/sra_small.csv` `assets/sra_large.csv` available for testing
*   Add vcf files generated from previous runs `--add_vcf_file` | `assets/precomputed/vcf` available for testing
*   Skip combined variant analysis (run reference prep and mapping) `skip_combined_analysis`
*   Skip samples capability `--skip_samples`
*   Skip samples file capability `--skip_samples_file`

### `Fixed`

*   
### `Dependencies`

*   Nextflow
*   nf-core
### `Deprecated`

*   GeneFlow support
*   QC pdf report
### `TODO`

*   Intermediate file cleanup and management
*   Update logo and metro-style workflow
