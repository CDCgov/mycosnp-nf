# nf-core/mycosnp: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of CDCgov/mycosnp-nf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

*   Support for phylogenetic tree generation
*   Skip samples capability `--skip_samples`
*   Skip samples file capability `--skip_samples_file`
*   Skip combined variant analysis (run reference prep and mapping) `skip_combined_analysis`
*   Skip reference generation `--ref_dir`
*   `--add_vcf_file`
### `Fixed`

*   
### `Dependencies`

*   Nextflow
*   nf-core
### `Deprecated`

*   QC pdf report

### `TODO`

*   Intermediate file cleanup and management