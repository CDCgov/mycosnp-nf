[![Open CDCgov/mycosnp-nf in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/CDCgov/mycosnp-nf)

Once the pod launches, it will present a VS-Code interface and comes with Nextflow, Conda and Docker pre-installed

* To run the pipeline with test data
```bash
nextflow run main.nf -profile docker,test
```

# ![CDCgov/mycosnp-nf](docs/images/nf-core-mycosnp_logo_light.png#gh-light-mode-only) ![CDCgov/mycosnp-nf](docs/images/nf-core-mycosnp_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/CDCgov/mycosnp-nf/workflows/nf-core%20CI/badge.svg)](https://github.com/CDCgov/mycosnp-nf/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/CDCgov/mycosnp-nf/workflows/nf-core%20linting/badge.svg)](https://github.com/CDCgov/mycosnp-nf/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/mycosnp/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23mycosnp-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/mycosnp)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction


**nf-core/mycosnp** is a bioinformatics best-practice analysis pipeline for MycoSNP is a portable workflow for performing whole genome sequencing analysis of fungal organisms, including Candida auris. This method prepares the reference, performs quality control, and calls variants using a reference. MycoSNP generates several output files that are compatible with downstream analytic tools, such as those for used for phylogenetic tree-building and gene variant annotations..

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/mycosnp/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

### Reference Preparation

> **Prepares a reference FASTA file for BWA alignment and GATK variant calling by masking repeats in the reference and generating the BWA index.**
* Genome repeat identification and masking (`nucmer`)
* BWA index generation (`bwa`)
* FAI and DICT file creation (`Picard`, `Samtools`)

### Sample QC and Processing

> **Prepares samples (paired-end FASTQ files) for GATK variant calling by aligning the samples to a BWA reference index and ensuring that the BAM files are correctly formatted. This step also provides different quality reports for sample evaluation.**

* Combine FASTQ file lanes if they were provided with multiple lanes.
* Filter unpaired reads from FASTQ files (`SeqKit`).
* Down sample FASTQ files to a desired coverage or sampling rate (`SeqTK`).
* Trim reads and assess quality (`FaQCs`).
* Generate a QC report by extracting data from FaQCs report data.
* Align FASTQ reads to a reference (`BWA`).
* Sort BAM files (`SAMTools`).
* Mark and remove duplicates in the BAM file (`Picard`).
* Clean the BAM file (`Picard "CleanSam"`).
* Fix mate information in the BAM file (`Picard "FixMateInformation"`).
* Add read groups to the BAM file (`Picard "AddOrReplaceReadGroups"`).
* Index the BAM file (`SAMTools`).
* [FastQC](#fastqc) - Filtered reads QC.
* Qualimap mapping quality report.
* [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline

### Variant calling and analysis

> **Calls variants and generates a multi-FASTA file and phylogeny.**

* Call variants (`GATK HaplotypeCaller`).
* Combine gVCF files from the HaplotypeCaller into a single VCF (`GATK CombineGVCFs`).
* Call genotypes using the (`GATK GenotypeGVCFs`).
* Filter the variants (`GATK VariantFiltration`) [default (but customizable) filter: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10'].
* Run a customized VCF filtering script (`Broad Institute`).
* Split the filtered VCF file by sample.
* Select only SNPs from the VCF files (`GATK SelectVariants`).
* Split the VCF file with SNPs by sample.
* Create a consensus sequence for each sample (`BCFTools`, `SeqTK`).
* Create a multi-fasta file from the VCF SNP positions using a custom script (`Broad`).
* Create phylogeny from multi-fasta file (`rapidNJ`, `FastTree2`, `RaxML`, `IQTree`)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run CDCgov/mycosnp-nf -profile test,YOURPROFILE
    ```

    Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

    > * The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` and are persistently observing issues downloading Singularity images directly due to timeout or network issues, then you can use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, you can use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```console
    nextflow run CDCgov/mycosnp-nf -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --fasta c_auris.fasta
    ```

## Documentation

The nf-core/mycosnp pipeline comes with documentation about the pipeline [usage](https://github.com/CDCgov/mycosnp-nf/blob/master/docs/usage.md), [parameters](https://github.com/CDCgov/mycosnp-nf/wiki/Parameter-Docs) and [output](https://github.com/CDCgov/mycosnp-nf/blob/master/docs/output.md).

## Credits

nf-core/mycosnp was originally written by CDC.

We thank the following people for their extensive assistance in the development of this pipeline:

* Michael Cipriano [@mjcipriano](https://github.com/mjcipriano)
* Sateesh Peri [@sateeshperi](https://github.com/sateeshperi)
* Hunter Seabolt [@hseabolt](https://github.com/hseabolt)
* Chris Sandlin [@cssandlin](https://github.com/cssandlin)
* Drewry Morris [@drewry](https://github.com/drewry)
* Lynn Dotrang [@leuthrasp](https://github.com/LeuThrAsp)
* Christopher Jossart [@cjjossart](https://github.com/cjjossart)
* Robert A. Petit III [@rpetit3](https://github.com/rpetit3)
## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  CDCgov/mycosnp-nf for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

# CDCgov GitHub Organization Open Source Project Template

**Template for clearance: This project serves as a template to aid projects in starting up and moving through clearance procedures. To start, create a new repository and implement the required [open practices](open_practices.md), train on and agree to adhere to the organization's [rules of behavior](rules_of_behavior.md), and [send a request through the create repo form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) using language from this template as a Guide.**

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 


## Related documents

* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)

  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
