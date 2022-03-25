# CDCgov/mycosnp-nf: Output

## Introduction

*   MycoSNP is a portable workflow for performing whole genome sequencing analysis of fungal organisms, including Candida auris.
*   This method prepares the reference, performs quality control, and calls variants using a reference.
*   MycoSNP generates several output files that are compatible with downstream analytic tools, such as those for used for phylogenetic tree-building and gene variant annotations.

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.


```bash
results/
├── combined
├── fastq
├── input
├── multiqc
├── pipeline_info
├── qc
├── reference
├── samples
└── stats
```

## Pipeline Overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [CDCgov/mycosnp: Output](#cdcgov-mycosnp-output)
	- [Introduction](#introduction) 
	- [Pipeline Overview](#pipeline-overview)
- [BWA Reference](#bwa-reference)
	- [Reference Preparation](#reference-preparation)
- [BWA Pre-process](#bwa-pre-process)
	- [Sample QC and Processing](#sample-qc-and-processing)
- [GATK Variants](#gatk-variants)
	- [Variant calling and analysis](#variant-calling-and-analysis)
- [Summary Files](#summary-files) 
	- [FastQC](#fastqc)
	- [QC report](#qc-report)
	- [MultiQC](#multiqc)
- [Pipeline Information](#pipeline-information)

## BWA Reference

### Reference Preparation

<details markdown="1">
<summary>Output files</summary>

* `reference/bwa/bwa`
    * `reference.amb`
    * `reference.ann`
    * `reference.bwt`
    * `reference.pac`
    * `reference.sa`
* `reference/dict`
	* `reference.dict`
* `reference/fai`
	* `reference.fa.fai`
* `reference/masked`
	* `reference.fa`    

</details>

> **Prepares a reference FASTA file for BWA alignment and GATK variant calling by masking repeats in the reference and generating the BWA index.**
* Genome repeat identification and masking (`nucmer`)
* BWA index generation (`bwa`)
* FAI and DICT file creation (`Picard`, `Samtools`)

## BWA Pre-process

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

**Important output files from this section:**

| File          | Path                                      |
| ---           | ---                                       |
| Trimmed Reads |  `(samples/<sampleID>/faqcs/*.fastq.gz)`  |
|  Bam files    |  `(samples/<sampleID>/finalbam)`          |
|  QC Report    |  `(stats/qc_report)`                      |
|  MultiQC      |  `(multiqc/)`                             |

## GATK Variants

### Variant calling and analysis

> **Calls variants, generates a multi-FASTA file, and creates phylogeny.**

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


**Important output files from this section:**

| File                                      | Path                                         |
| ---                                       | ---                                          |
| Individual VCF Files from HaplotypeCaller |  `(samples/variant_calling/haplotypecaller)` |
|  Filtered selected variants combined      |  `(combined/finalfiltered)`                  |
|  Filtered selected variants individual    |  `(combined/splitvcf)`                       |
|  Individual consensus fasta files         |  `(combined/consensus)`                      |
|  Selected SNP fasta file                  |  `(combined/vcf-to-fasta)`                   |
|  Phylogeny files                          |  `(combined/phylogeny/)`                     |


## Summary Files

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### QC Report
The QC report values are generated from FAQCS text file outputs. The following is an example table:
| Sample Name | # Reads Before Trimming | GC Before Trimming | Average Phred Before Trimming | Coverage Before Trimming | # Reads After Trimming | # Paired Reads After Trimming | # Unpaired Reads After Trimming | GC After Trimming | Average Phred After Trimming | Coverage After Trimming |
|-------------|-------------------------|--------------------|-------------------------------|--------------------------|------------------------|-------------------------------|---------------------------------|-------------------|------------------------------|-------------------------|
| ERR2172265  | 367402                  | 52.44%             | 34.64                         | 16.578758543095503       | 367396                 | 367390 (100.00 %)             | 6 (0.00 %)                      | 52.45%            | 34.64                        | 16.578487797287757      |

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
