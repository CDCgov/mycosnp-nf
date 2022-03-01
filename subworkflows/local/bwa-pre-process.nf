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
include { SAMTOOLS_VIEW as PICARDDUPTOCLEANSAM}                 from '../../modules/nf-core/modules/samtools/view/main'
include { PICARD_FIXMATEINFORMATION }     from '../../modules/nf-core/modules/picard/fixmateinformation/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/nf-core/modules/picard/addorreplacereadgroups/main'
include { SAMTOOLS_INDEX }                from '../../modules/nf-core/modules/samtools/index/main'
include { FASTQC }                        from '../../modules/nf-core/modules/fastqc/main'
include { MULTIQC }                       from '../../modules/nf-core/modules/multiqc/main'
include { QUALIMAP_BAMQC }                from '../../modules/nf-core/modules/qualimap/bamqc/main'
include { DOWNSAMPLE_RATE }               from '../../modules/local/downsample_rate.nf'



workflow BWA_PREPROCESS {

    take:
    reference //channel: [ tuple reference_fasta, samtools_faidx, bwa_index ]
    reads // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions           = Channel.empty()
    ch_alignment          = Channel.empty()
    ch_alignment_index    = Channel.empty()
    ch_alignment_combined = Channel.empty()
    

    // Going to skip this one for now and add multilane support within the read_csv
    // CONCAT_FASTQ_LANES()
    
    SEQKIT_PAIR(reads)
	DOWNSAMPLE_RATE(SEQKIT_PAIR.out.reads, reference[0], params.coverage)
	DOWNSAMPLE_RATE.out.downsampled_rate.view()
    SEQTK_SAMPLE(SEQKIT_PAIR.out.reads, DOWNSAMPLE_RATE.out.number_to_sample)
    FAQCS(SEQTK_SAMPLE.out.reads)
    // QC() // Local qc report
    BWA_MEM(FAQCS.out.reads, reference[2], true) // This already sorts the bam file
    PICARD_MARKDUPLICATES(BWA_MEM.out.bam)
    PICARD_CLEANSAM(PICARD_MARKDUPLICATES.out.bam)
    PICARD_FIXMATEINFORMATION(PICARD_CLEANSAM.out.bam)
    PICARD_ADDORREPLACEREADGROUPS(PICARD_FIXMATEINFORMATION.out.bam)
    SAMTOOLS_INDEX(PICARD_ADDORREPLACEREADGROUPS.out.bam)
    FASTQC(reads)
    // QUALIMAP_BAMQC()
    // MULTIQC()

    ch_alignment_combined = PICARD_ADDORREPLACEREADGROUPS.out.bam.join(SAMTOOLS_INDEX.out.bai)
    
  
    ch_versions            = ch_versions.mix(  SEQKIT_PAIR.out.versions, 
                                               SEQTK_SAMPLE.out.versions, 
                                               FAQCS.out.versions,
                                               BWA_MEM.out.versions,
                                               PICARD_MARKDUPLICATES.out.versions,
                                               PICARD_CLEANSAM.out.versions,
                                               PICARD_FIXMATEINFORMATION.out.versions,
                                               PICARD_ADDORREPLACEREADGROUPS.out.versions,
                                               SAMTOOLS_INDEX.out.versions,
                                               FASTQC.out.versions
                                               //,
                                               //QUALIMAP_BAMQC.out.versions,
                                               //MULTIQC.out.versions
                                            )
    ch_alignment          = PICARD_ADDORREPLACEREADGROUPS.out.bam
    ch_alignment_index    = SAMTOOLS_INDEX.out.bai


    emit:
    alignment          = ch_alignment          // channel: [ val(meta), [ bam ] ]
    alignment_index    = ch_alignment_index    // channel: [ val(meta), [ bai ] ]
    alignment_combined = ch_alignment_combined // channel: [ val(meta), [ vcf ] ]
    versions              = ch_versions        // channel: [ ch_versions ]
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
