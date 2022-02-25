/*
========================================================================================
    BWA Pre-Process Sub-Workflow
========================================================================================
*/

include { SEQKIT_PAIR }                   from '../../modules/nf-core/modules/seqkit/pair/main'
include { SEQTK_SAMPLE }                  from '../../modules/nf-core/modules/seqtk/sample/main'
include { FAQCS }                         from '../../modules/nf-core/modules/faqcs/main'
// TODO: QC report local module
include { BWA_INDEX }                     from '../../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM }                       from '../../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_SORT }                 from '../../modules/nf-core/modules/samtools/sort/main'
include { PICARD_MARKDUPLICATES }         from '../../modules/nf-core/modules/picard/markduplicates/main'
include { PICARD_CLEANSAM }               from '../../modules/nf-core/modules/picard/cleansam/main'
include { PICARD_FIXMATEINFORMATION }     from '../../modules/nf-core/modules/picard/fixmateinformation/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/nf-core/modules/picard/addorreplacereadgroups/main'
include { SAMTOOLS_INDEX }                from '../../modules/nf-core/modules/samtools/index/main'
include { FASTQC }                        from '../../modules/nf-core/modules/fastqc/main'
include { MULTIQC }                       from '../../modules/nf-core/modules/multiqc/main'
include { QUALIMAP_BAMQC }                from '../../modules/nf-core/modules/qualimap/bamqc/main'



workflow BWA_PREPROCESS {

    take:   
    reference //channel: tuple reference_fasta, samtools_faidx, bwa_index
    reads // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()


    // CONCAT_FASTQ_LANES()
    // SEQKIT_PAIR()
    // SEQTK_SAMPLE()
    // FAQCS(DOWNSAMPLE.out)
    // QC()
    BWA_MEM(reads, reference[2], true)
    // BWA_SORT()
    // MARK_DUPLICATES()
    // CLEAN_SAM()
    // FIXMATEINFORMATION()
    // ADDORREPLACEGROUPS()
    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    FASTQC(reads)
    // QUALIMAP()
    // MULTIQC()
    ch_combined = Channel.empty()
    BWA_MEM.out.bam.combine(SAMTOOLS_INDEX.out.bai).map{meta1, bam, meta2, bai -> [meta1, bam, bai] }.set{ch_combined}

    ch_versions = ch_versions.mix(  BWA_MEM.out.versions, 
                                    SAMTOOLS_INDEX.out.versions, 
                                    FASTQC.out.versions
                                 )

    emit:
    //alignment = ADDORREPLACEGROUPS.out
    alignment = BWA_MEM.out.bam
    alignment_index = SAMTOOLS_INDEX.out.bai
    alignment_combined = ch_combined
    versions = ch_versions // channel: [ versions.yml ]

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
