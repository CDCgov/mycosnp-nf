# CDCgov/mycosnp-nf: Output

## Introduction
This document describes the output produced by both the Pre-MycoSNP and main MycoSNP workflows.

The directories listed below will be created in the results directory after either workflow has finished. All paths are relative to the top-level results directory.

```bash
results/
├── combined
├── fastq
├── input
├── multiqc
├── pipeline_info
├── reference
├── samples
├── snpeff
└── stats
```

## Navigation
- [Pre-MycoSNP workflow](#pre-mycosnp-workflow-overview)
- [Main MycoSNP workflow](#main-mycosnp-workflow-overview)

## Pre-MycoSNP workflow: Overview
> [!NOTE]
> Pre-MycoSNP's taxonomic classification and subtyping steps (for _C. auris_ clade typing) have been validated for common fungal pathogens. See the [report](https://github.com/CDCgov/mycosnp-nf/wiki/Validation-study:-Pre%E2%80%90MycoSNP-taxonomic-classification-and-subtyping) on the Wiki for more.

- [Quality control](#quality-control)
    - [Sequencing read quality metrics (FastQC)](#sequencing-read-quality-metrics-fastqc)
    - [Read trimming/filtering (FaQCS)](#read-trimmingfiltering-faqcs)
- [De novo assembly (Shovill)](#de-novo-assembly-shovill)
- [Taxonomic classification (GAMBIT)](#taxonomic-classification-gambit)
- [Subtyping (sourmash)](#subtyping-sourmash)
- [Summary Files](#summary-files)
  - [Pre-MycoSNP summary report](#pre-mycosnp-summary-report)
  - [MultiQC](#multiqc)
  - [Pipeline information](#pipeline-information)

## Quality control
### Sequencing read quality metrics (FastQC)
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides general quality metrics about sequenced reads, e.g Q scores, GC content, adapter contamination.

Output files:
- `samples/<sample_id>/fastqc_raw/`

### Read trimming/filtering (FaQCs)
[FaQCs](https://github.com/LANL-Bioinformatics/FaQCs) trims reads and assesses quality.

Output files:
- `samples/<sample_id>/faqcs/`

## De novo assembly (Shovill)
- [Shovill](https://github.com/tseemann/shovill) performs de novo assembly using an assembler of your choice ([SKESA](https://github.com/ncbi/SKESA) [default in the Pre-MycoSNP workflow], [SPAdes](https://github.com/ablab/spades), [Megahit](https://github.com/voutcn/megahit), or [Velvet](https://github.com/dzerbino/velvet​)), handling downsampling and read trimming as part of the pipeline.
- By default, Pre-MycoSNP sets Shovill to downsample to 70X depth.
- To save time, and because high-quality assemblies aren't necessary for Pre-MycoSNP's sketch-based taxonomic classification and clade typing, Pre-MycoSNP skips Shovill's correction step where reads are mapped back to contigs to correct minor assembly errors (by setting the `--nocorr` flag).

Output files:
- `samples/<sample_id>/assembly/<sample_id>.fa`: Assembly in FASTA format

## Taxonomic classification (GAMBIT)
- Pre-MycoSNP uses [GAMBIT](https://github.com/jlumpe/gambit) to perform taxonomic classification from de novo assemblies.
- GAMBIT makes compressed sketches of the k-mer content of an assembly and compares it to sketches in a database.
- It classifies samples to the species level if possible. If not possible, it attempts to classify to the genus level. If not able to classify to the genus or species levels, it will not make a prediction, but will still report the closest match in the database.
- GAMBIT reports distances between two sketches with the Jaccard distance (1 - [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index)).
- Pre-MycoSNP uses GAMBIT's fungal database v1.0.0. See the [list of taxa included in the database](/assets/gambit_db/gambit-1.0.0-20241213-taxa-list.txt). See [GAMBIT'S documentation](https://theiagen.notion.site/GAMBIT-7c1376b861d0486abfbc316480046bdc#3f6610c81fbb4812b745234441514e12) for more details about the database.

Output files:
- `samples/<sample_id>/taxonomy/<sample_id>_gambit.csv`: GAMBIT output

## Subtyping (sourmash)
- After the GAMBIT step, Pre-MycoSNP uses [sourmash](https://github.com/sourmash-bio/sourmash) to perform subtyping for taxa that are contained in [`assets/sourmash_db/sourmash_taxa.csv`](/assets/sourmash_db/sourmash_taxa.csv) and have an associated sourmash signature file in [`assets/sourmash_db/signatures/`](/assets/sourmash_db/signatures/).
- For _Candida auris_, clade prediction (Clades I-VI) is performed with [`assets/sourmash_db/signatures/candida_auris_clades.sig`](/assets/sourmash_db/signatures/candida_auris_clades.sig).
- Pre-MycoSNP estimates average nucleotide identity (ANI) to the closest match in the signature file using the `--containment` and `--ani` options in [`sourmash compare`](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-compare-compare-many-signatures).

Output files:
- `samples/<sample_id>/subtype/<sample_id>_subtype.csv`: Subtyping output. Empty if subtyping not performed.

## Summary Files
### Pre-MycoSNP summary report
Output files: `combined/pre-mycosnp_summary/pre-mycosnp-summary.csv`

#### Header descriptions and comparison to main MycoSNP workflow's [QC report](#qc-report-statsqc_reportqc_reporttxt):
> [!NOTE]
> If a _C. auris_ sample has a `Subtype_ANI` less than 99.7, the `Subtype_Closest_Match` field is filled with "ANI is less than the established Candida auris clade separation threshold of 99.7." This is based on reliable clade separation thresholds described in the [validation report](https://github.com/CDCgov/mycosnp-nf/wiki/Validation-study:-Pre%E2%80%90MycoSNP-taxonomic-classification-and-subtyping) on the Wiki. _C. auris_ clade predictions with ANI less than 99.7 should be interpreted with caution.

| Column                           | Description                                                                                 | Matching column in main MycoSNP workflow's QC report |
| -------------------------------- | ------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| Sample                           | sample ID                                                                                   | `Sample Name`                                        |
| Predicted_Rank                   | Output from GAMBIT: "species", "genus", or "no prediction"                                  |                                                      |
| Predicted_Taxon                  | GAMBIT prediction                                                                           |                                                      |
| Subtype_Closest_Match            | sourmash closest match (see "Note" below)                                                   |                                                      |
| Subtype_ANI                      | ANI to closest sourmash match (see "Note" below)                                            |                                                      |
| Closest_GAMBIT_Entry_Description | Closest entry in the GAMBIT fungal database                                                 |                                                      |
| Closest_GAMBIT_Entry_Distance    | Jaccard distance (`1 - Jaccard index`) to closest entry                                     |                                                      |
| Trimmed_Reads                    | Number of reads after trimming with FaQCs                                                   | `Reads After Trimming`                               |
| Avg_Read_Quality                 | Average Q score of reads after trimming with FaQCs                                          | `Average Q Score After Trimming`                     |
| Sample_Assembly_Length           | Length of the Shovill assembly                                                              |                                                      |
| Sample_Assembly_GC               | GC percentage of Shovill assembly                                                           |                                                      |
| Reference_Genome_Length          | Length of the reference genome of the closest GAMBIT entry (if a prediction is made)        |                                                      |
| Avg_Depth_Coverage               | (Total trimmed bases / `Reference_Genome_Length`)                                           |                                                      |
| Reference_GC                     | GC percentage of the reference genome of the closest GAMBIT entry (if a prediction is made) |                                                      |

### MultiQC
- Outputs follow the same [pattern described for the main MycoSNP workflow](#multiqc-1).
- Output files: `multiqc/pre-mycosnp/`

### Pipeline information
Outputs follow the same [pattern described for the main MycoSNP workflow](#pipeline-information-1).

## Main MycoSNP workflow: Overview
- [BWA Reference](#bwa-reference)
    - [Reference Preparation](#reference-preparation)
- [BWA Pre-process](#bwa-pre-process)
  - [Sample QC and Processing](#sample-qc-and-processing)
- [GATK Variants](#gatk-variants)
  - [Variant calling and analysis](#variant-calling-and-analysis)
- [Variant Annotation](#variant-annotation)
  - [snpEff analysis](#snpeff-analysis)
- [Summary Files](#summary-files-1)
  - [FastQC](#fastqc)
  - [QC Report](#qc-report-statsqc_reportqc_reporttxt)
  - [MultiQC](#multiqc-1)
  - [Pipeline information](#pipeline-information-1)

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

| File                                               | Path                                      |
| ---                                                | ---                                       |
| Trimmed Reads                                      |  `(samples/<sample_id>/faqcs/*.fastq.gz)` |
| Bam files                                          |  `(samples/<sample_id>/finalbam)`         |
| [QC Report](#qc-report-statsqc_reportqc_reporttxt) |  `(stats/qc_report/qc_report.txt)`        |
| [MultiQC](#multiqc-1)                              |  `(multiqc/)`                             |

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
* Create a multi-fasta file from the VCF SNP positions using a custom script (`Broad`).
* Create distance matrix file using muti-fasta file (`snp-dists`).
* Create phylogeny from multi-fasta file (`rapidNJ`, `FastTree2`, `quicksnp`,`RaxML(optional)`, `IQTree(optional)`)


**Important output files from this section:**

| File                                      | Path                                                    |
| ---                                       | ---                                                     |
| Individual VCF Files from HaplotypeCaller |  `(samples/<sample_id>variant_calling/haplotypecaller)` |
|  Filtered selected variants combined      |  `(combined/finalfiltered)`                             |
|  Filtered selected variants individual    |  `(combined/splitvcf)`                                  |
|  Selected SNP fasta file                  |  `(combined/vcf-to-fasta)`                              |
|  Phylogeny files                          |  `(combined/phylogeny/)`                                |


## Variant Annotation 
### snpEff analysis 
> **Annotate variants within the FKS1 gene using the filtered vcf file and provide output in the form of a report (currently for _C. auris_ B11205 reference only)**

* snpEff annotation (`snpeff`)
* [`snpeffr`](https://github.com/CDCgov/snpeffr) creates a report using the snpeff annotated vcf file. Non-synonymous variants in FKS1 hotspot regions are included in the report (snpeff/combined_cauris_refB11205_fks1.csv)."

**Important output files from this section:**

| File                                           | Path                                          |
| ---                                            | ---                                           |
| Annotated vcf file from SnpEff	               | `(snpeff/combined.snpeff.vcf.gz)`             |
| Customized report of variants in FKS1 hotspots | `(snpeff/combined_cauris_refB11205_fks1.csv)` |
> [!NOTE]
> The "mutation" column in `snpeff/combined_cauris_refB11205_fks1.csv` can contain a value of "undetermined". This occurs when GATK identifies a variant in the individual VCF file for an isolate, but then the identified variant is lost in the subsequent step when GATK creates a merged VCF file that includes all isolates.

## Summary Files

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> [!NOTE]
> The FastQC plots displayed in the MultiQC report show _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### QC Report (`stats/qc_report/qc_report.txt`)
The QC report values are generated from FAQCS text file outputs and Qualimap result files. The following is an example table:
| Sample Name | Reads Before Trimming | GC Before Trimming | Average Q Score Before Trimming | Reference Length Coverage Before Trimming | Reads After Trimming | Paired Reads After Trimming | Unpaired Reads After Trimming | GC After Trimming | Average Q Score After Trimming | Reference Length Coverage After Trimming | Mean Coverage Depth[^1] | Reads Mapped     | Genome Fraction at 10X[^2] |
|-------------|-----------------------|--------------------|---------------------------------|-------------------------------------------|----------------------|-----------------------------|-------------------------------|-------------------|--------------------------------|------------------------------------------|---------------------|------------------|----------------------------|
| SRR10461157 | 3987618               | 44.73%             | 35.03                           | 80.68                                     | 3987537 (100.00 %)   | 3987456 (100.00 %)          | 81 (0.00 %)                   | 44.73%            | 35.03                          | 80.68                                    | 74.66               | 3789190 (95.31%) | 97.70%                     |

| QC Metric                                  | Source   |
|--------------------------------------------|----------|
| Reads Before Trimming                      | FAQCS    |
| GC Before Trimming                         | FAQCS    |
| Average Q Score Before Trimming            | FAQCS    |
| Reference Length Coverage Before Trimming  | FAQCS    |
| Reads After Trimming                       | FAQCS    |
| Paired Reads After Trimming                | FAQCS    |
| Unpaired Reads After Trimming              | FAQCS    |
| GC After Trimming                          | FAQCS    |
| Average Q Score After Trimming             | FAQCS    |
| Reference Length Coverage After Trimming   | FAQCS    |
| Mean Coverage Depth[^1]                    | Qualimap |
| Reads Mapped                               | Qualimap |
| Genome Fraction at 10X[^2]                 | Qualimap |

[^1]: "Mean Coverage Depth" is a post-alignment metric. Qualimap divides the reference genome into 400 (default) windows. For each window, Qualimap counts the total number of bases that aligned with the window, divided by the number of bases in the window. The Mean Coverage Depth is the average across all of these windows.
[^2]: "Genome fraction" refers to the percentage of the reference genome that is covered at a specific depth by the sample. This metric will output the genome fraction at the minimum depth (`--min_depth`) parameter used to call bases (e.g. "Genome Fraction at 10X").

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/mycosnp/`
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

Nextflow provides excellent functionality for generating various [reports](https://www.nextflow.io/docs/latest/reports.html) relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
