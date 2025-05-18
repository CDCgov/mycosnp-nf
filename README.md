# 🍄🧬 MycoSNP: Whole Genome Sequencing Analysis of Fungal Isolates

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction
**MycoSNP** is a bioinformatics pipeline for performing whole genome sequencing analysis of fungal organisms (e.g. _Candida auris_) from Illumina paired-end reads. It is built with the [nf-core](https://nf-co.re/) template.

This repository contains two workflows that are run independently:
- **Pre-MycoSNP workflow**: A first-pass workflow for quick answers
  - Fungal taxonomic classification and _Candida auris_ clade typing, using de novo assemblies
- **Main MycoSNP workflow (default workflow)**:
  - Reference-based SNP calling
  - Tree building
  - Identification of antifungal-resistance mutations

## 📚 Full Documentation
- [**Usage**](docs/usage.md): An overview of how MycoSNP works, how to run it and a description of all of the different command-line flags.
- [**Parameters**](docs/params.md): Options, flags, and inputs.
- [**Output**](docs/output.md): An overview of the different results produced by MycoSNP and how to interpret them.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_. More info on using containers with Nextflow [here](https://www.nextflow.io/docs/latest/container.html#container-page).

> [!TIP]
> Using Apptainer/Singularity with Nextflow version >=23 can result in failures in Linux server environments due to peculiarities with container directory mounting. If you are experiencing `No such file or directory` errors, try running with an earlier version of Nextflow (we've had success with 22.10.6).
3. Test the main MycoSNP workflow on pre-defined minimal test samples with a single command:
    ```console
    nextflow run CDCgov/mycosnp-nf -profile test,YOURPROFILE
    ```
> [!NOTE]
> The samples for the test run are bacterial (_N. gonorrhoeae_), not fungal. This is intentional so the test finishes in a few minutes (as opposed to longer for fungal samples with much larger genomes).

  > [!TIP]
  > Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.
  > * The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
  > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
  > * If you are using `singularity` and are persistently observing issues downloading Singularity images directly due to timeout or network issues, then you can use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, you can use the [`nf-core download`](https://nf-co.re/docs/nf-core-tools/pipelines/download) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/container.html#singularity) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
  > * If you are using `conda` (not recommended), it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!
> [!NOTE]
> The `--workflow` option specifies which workflow to run (Pre-MycoSNP workflow or main MycoSNP workflow). By default (when no `workflow` is specified), the main MycoSNP workflow is executed. See summaries of each workflow in the next section, or more detailed documentation in the [Full Documentation](docs/README.md).
  * [Pre-MycoSNP workflow](#pre-mycosnp-workflow-summary):
    ```console
    nextflow run CDCgov/mycosnp-nf --workflow PRE_MYCOSNP -profile <docker/singularity/other/institute> --input samplesheet.csv
    ```
  * [Main MycoSNP workflow](#main-mycosnp-workflow-default-workflow-summary) (default workflow):
    ```console
    nextflow run CDCgov/mycosnp-nf -profile <docker/singularity/other/institute> --input samplesheet.csv --fasta reference_genome.fasta
    ```
5. It is advisable to delete large temporary or log files after the successful completion of the run. It takes a lot of space and may cause issues in future runs.

## Pre-MycoSNP Workflow: Summary
1. Sequencing reads quality metrics ([`FastQC`](https://github.com/s-andrews/FastQC))
1. Remove unpaired reads ([`seqkit pair`](https://github.com/shenwei356/seqkit))​
1. Trim/filter reads ([`FaQCs`](https://github.com/LANL-Bioinformatics/FaQCs))​
1. De novo assembly ([`Shovill`](https://github.com/tseemann/shovill))
    - Downsampling​ (default 70X)
    - Assembly (with [`SKESA`](https://github.com/ncbi/SKESA) [default], [`SPAdes`](https://github.com/ablab/spades), [`Megahit`](https://github.com/voutcn/megahit), or [`Velvet`](https://github.com/dzerbino/velvet​)​)
1. Taxonomic classification ([`GAMBIT`](https://github.com/jlumpe/gambit​))
    - Classifies isolate to genus/species level, if possible
    - Uses GAMBIT's fungal database v1.0.0. See the [list of taxa included in the database](assets/gambit_db/gambit-1.0.0-20241213-taxa-list.txt). See [GAMBIT'S documentation](https://theiagen.notion.site/GAMBIT-7c1376b861d0486abfbc316480046bdc#3f6610c81fbb4812b745234441514e12) for more details about the database.
1. Subtyping ([`sourmash`](https://github.com/sourmash-bio/sourmash))
    - Compares sourmash sketch of sample against sourmash signature file provided in [`assets/sourmash_db/`](assets/sourmash_db/).
    - By default, for _C. auris_, Pre-MycoSNP performs clade typing (Clades I-VI).
1. [Summary report](docs/output.md#pre-mycosnp-summary-report) containing genus/species classification, subtype, and read and assembly quality metrics.

## Main MycoSNP Workflow (Default Workflow): Summary

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
* FastQC - Filtered reads QC.
* Qualimap mapping quality report.
* MultiQC - Aggregate report describing results and QC from the whole pipeline

### Variant calling and analysis

> **Calls variants and generates a multi-FASTA file and phylogeny.**

* Call variants (`GATK HaplotypeCaller`).
* Combine gVCF files from the HaplotypeCaller into a single VCF (`GATK CombineGVCFs`).
* Call genotypes using the (`GATK GenotypeGVCFs`).
* Filter the variants (`GATK VariantFiltration`) [default (but customizable) filter: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10'].
* Run a customized VCF filtering script (Broad Institute).
* Split the filtered VCF file by sample.
* Select only SNPs from the VCF files (`GATK SelectVariants`).
* Split the VCF file with SNPs by sample.
* Create a multi-fasta file from the VCF SNP positions using a custom script (Broad Institute).
* Create a distance matrix file using multi-fasta file(`SNPdists`).
* Create phylogeny from multi-fasta file (`rapidNJ`, `FastTree2`, `quicksnp`, `RaxML(optional)`, `IQTree(optional)`)

### Variant annotation analysis (currently available for *C. auris* B11205 genome only)

* annotated VCF file (`snpEff`)
* [`snpeffr`](https://github.com/CDCgov/snpeffr) report. Non-synonymous variants in FKS1 hotspot regions are included in the report.

## Pre-configured Nextflow development environment using Gitpod

>[![Open CDCgov/mycosnp-nf in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/CDCgov/mycosnp-nf)
>>Once the pod launches, it will present a VS-Code interface and comes with Nextflow, Conda and Docker pre-installed

## Credits

nf-core/mycosnp was originally developed at CDC (with contributions from many others in the public health bioinformatics community), and is currently maintained by CDC's Mycotic Diseases Branch.

We thank the following people (alphabetical order by last name) for their code contributions or assistance in the development of this pipeline:

* John Arnn [@jwarnn](https://github.com/jwarnn)
* Ujwal Bagal [@urbagal](https://github.com/urbagal)
* Arun Boddapati [@arunbodd](https://github.com/arunbodd)
* Michael Cipriano [@mjcipriano](https://github.com/mjcipriano)
* Lynn Dotrang [@leuthrasp](https://github.com/LeuThrAsp)
* Jared Johnson [@DOH-JDJ0303](https://github.com/DOH-JDJ0303)
* Christopher Jossart [@cjjossart](https://github.com/cjjossart)
* Elizabeth Misas [@AspTryGlu](https://github.com/AspTryGlu)
* Drewry Morris [@drewry](https://github.com/drewry)
* Zack Mudge [@zmudge3](https://github.com/zmudge3)
* Harshil Patel [@drpatelh](https://github.com/drpatelh)
* Sateesh Peri [@sateeshperi](https://github.com/sateeshperi)
* Robert A. Petit III [@rpetit3](https://github.com/rpetit3)
* Malavika Rajeev [@mrajeev08](https://github.com/mrajeev08)
* Charlotte Royer [@royercj](https://github.com/royercj)
* Chris Sandlin [@cssandlin](https://github.com/cssandlin)
* Hunter Seabolt [@hseabolt](https://github.com/hseabolt)

> Special thanks to [**StaPH-B**](https://staphb.org/) for open-source collaborations and discussions.

> Special thanks to CDC's Office of Advanced Molecular Detection (OAMD) and the Scientific Computing and Bioinformatics Support (SciComp) team for supporting development and computing infrastructure.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `MycoSNP` and `nf-core` publications as follows:

> Bagal UR, Phan J, Welsh RM, Misas E, Wagner D, Gade L, Litvintseva AP, Cuomo CA, Chow NA. 
> 
> MycoSNP: A Portable Workflow for Performing Whole-Genome Sequencing Analysis of Candida auris. 
> 
> _Methods Mol Biol._ 2022; 2517:215-228. doi: [10.1007/978-1-0716-2417-3_17](https://pubmed.ncbi.nlm.nih.gov/35674957)

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

# CDCgov GitHub Organization Open Source Project

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Related documents

* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)
  
## Public Domain
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.
## License
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

## Privacy
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Information
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
