/*
========================================================================================
    BWA Pre-Process Sub-Workflow
========================================================================================
*/

include { SEQKIT_PAIR }              from '../../modules/nf-core/modules/seqkit/pair/'
include { SEQTK_SAMPLE }              from '../../modules/nf-core/modules/seqtk/sample/'
include { FAQCS }                     from '../../modules/nf-core/modules/faqcs/'
// TODO: QC report local module
include { BWA_INDEX }                 from '../../modules/nf-core/modules/bwa/index/'
include { BWA_MEM }                   from '../../modules/nf-core/modules/bwa/mem/'
include { SAMTOOLS_SORT }             from '../../modules/nf-core/modules/samtools/sort/'
include { PICARD_MARKDUPLICATES }     from '../../modules/nf-core/modules/picard/markduplicates/'
include { PICARD_CLEANSAM }           from '../../modules/nf-core/modules/picard/cleansam/'
include { PICARD_FIXMATEINFORMATION } from '../../modules/nf-core/modules/picard/fixmateinformation/'
include { PICARD_ADDORREPLACEGROUPS } from '../../modules/nf-core/modules/picard/addorreplacegroups/'
include { SAMTOOLS_INDEX }            from '../../modules/nf-core/modules/samtools/index/'
include { FASTQC }                    from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC }                   from '../modules/nf-core/modules/multiqc/main'
include { QUALIMAP_BAMQC }            from '../modules/nf-core/modules/qualimap/bamqc/'



workflow BWA_PREPROCESS {

    take:
    tuple reference_fasta, samtools_faidx, bwa_index
    tuple meta, fastq

    main:
    // CONCAT_FASTQ_LANES()
    // FASTQ_PAIR()
    DOWNSAMPLE()
    FAQCS(DOWNSAMPLE.out)
    QC()
    BWA_ALIGN()
    BWA_SORT()
    MARK_DUPLICATES()
    CLEAN_SAM()
    FIXMATEINFORMATION()
    ADDORREPLACEGROUPS()
    BAM_INDEX()
    FASTQC()
    QUALIMAP()
    MULTIQC()


    emit:
    metaout = meta
    alignment = ADDORREPLACEGROUPS.out
    alignment_index = BAM_INDEX.out
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}

/*

Workflow Process:

1. Combine FASTQ file lanes if they were provided with multiple lanes.
2. Filter unpaired reads from FASTQ files using SeqKit 0.16.0.
seqkit pair -1 R1_FILE -2 R2_FILE
3. Down sample FASTQ files to a desired coverage or sampling rate using SeqTK 1.3.
seqt
4. Trim reads and assess quality using FaQCs 2.10.

5. Generate a QC report by extracting data from FaQCs PDFs.
custom linux stuff
6. Align FASTQ reads to a reference using BWA 0.7.17.
R1/R2 from Step 5
bwa mem -t THREADS REFERENCE_INDEX R1 R2
7. Sort BAM files using SAMTools 1.10.
8. Mark and remove duplicates in the BAM file using Picard 2.22.7.
9. Clean the BAM file using Picard 2.22.7 "CleanSam".
10. Fix mate information in the BAM file using Picard 2.22.7 "FixMateInformation".


11. Add read groups to the BAM file using Picard 2.22.7 "AddOrReplaceReadGroups".


12. Index the BAM file using SAMTools 1.10.
    cp ${INPUT_FULL} ${OUTPUT_FULL}/${OUTPUT_BASE}.bam
    samtools index ${OUTPUT} ${OUTPUT_BASE}/${OUTPUT_BASE}.bam


*/
