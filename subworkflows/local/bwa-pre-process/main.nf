/*
========================================================================================
    BWA Pre-Process Sub-Workflow
========================================================================================
*/

include { SEQKIT_PAIR                          } from '../../../modules/nf-core/seqkit/pair/main'
include { SEQTK_SAMPLE                         } from '../../../modules/nf-core/seqtk/sample/main'
include { FAQCS                                } from '../../../modules/local/faqcs/main'
include { BWA_INDEX                            } from '../../../modules/nf-core/bwa/index/main'   //not used here
include { BWA_MEM                              } from '../../../modules/local/bwa/mem/main'
include { SAMTOOLS_SORT                        } from '../../../modules/nf-core/samtools/sort/main' //not used here
include { PICARD_MARKDUPLICATES                } from '../../../modules/local/picard/markduplicates/main'
include { PICARD_CLEANSAM                      } from '../../../modules/nf-core/picard/cleansam/main'
include { SAMTOOLS_VIEW as PICARDDUPTOCLEANSAM } from '../../../modules/nf-core/samtools/view/main' //not used here
include { PICARD_FIXMATEINFORMATION            } from '../../../modules/nf-core/picard/fixmateinformation/main'
include { PICARD_ADDORREPLACEREADGROUPS        } from '../../../modules/local/picard/addorreplacereadgroups/main'
include { SAMTOOLS_INDEX                       } from '../../../modules/nf-core/samtools/index/main'
include { FASTQC as FASTQC_POST                } from '../../../modules/nf-core/fastqc/main'
include { QUALIMAP_BAMQC                       } from '../../../modules/nf-core/qualimap/bamqc/main'
include { DOWNSAMPLE_RATE                      } from '../../../modules/local/downsample_rate/main'
include { SAMTOOLS_STATS                       } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS                    } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT                    } from '../../../modules/nf-core/samtools/flagstat/main'
include { QC_REPORT                            } from '../../../modules/local/qc_report/main'



workflow BWA_PREPROCESS {

    take:
    ref_fasta // channel: tuple val(meta), path(fasta)
    ref_fai   // channel: tuple val(meta), path(fai)
    ref_bwa   // channel: tuple val(meta), path(bwa)
    reads     // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions           = Channel.empty()
    ch_alignment          = Channel.empty()
    ch_alignment_index    = Channel.empty()
    ch_alignment_combined = Channel.empty()
    ch_seq_samplerate     = Channel.empty()


    SEQKIT_PAIR(reads)

    if (params.coverage == 0) {
        FAQCS(SEQKIT_PAIR.out.reads)
        }
    else {
        // DOWNSAMPLE_RATE expects a plain path for reference_fasta
        DOWNSAMPLE_RATE(SEQKIT_PAIR.out.reads, ref_fasta.map{ meta, fa -> fa }, params.coverage)

        ch_seq_samplerate = SEQKIT_PAIR.out.reads.join(
            DOWNSAMPLE_RATE.out.downsampled_rate.map{
                meta, sr, snr -> [ meta, snr]
                }
            )

        // breakup file list into individual rows for processing
        ch_seq_samplerate
            .flatMap { meta, reads_list, size ->
                def (fastq1, fastq2) = reads_list
                [
                    tuple(meta, fastq1 ,size),
                    tuple(meta, fastq2, size)
                ]
            }
            .set { ch_seq_samplerate }

        SEQTK_SAMPLE( ch_seq_samplerate )

        // convert R1/R2 mapped reads back into single sorted channel rows
        ch_seq_sample_combine = SEQTK_SAMPLE.out.reads
            .groupTuple( by: 0, size: 2, sort: true )

        FAQCS( ch_seq_sample_combine )
        }

    // Group trimmed reads by sample so BWA_MEM receives [meta, [R1,R2]]
    // Pass only reads + index path + sort flag as required by local BWA_MEM (3 inputs)
    BWA_MEM(
        FAQCS.out.reads,
        ref_bwa.map{ meta, idx -> idx },
        true
    )

    PICARD_MARKDUPLICATES( BWA_MEM.out.bam )

    PICARD_CLEANSAM(PICARD_MARKDUPLICATES.out.bam )

    PICARD_FIXMATEINFORMATION( PICARD_CLEANSAM.out.bam )

    PICARD_ADDORREPLACEREADGROUPS( PICARD_FIXMATEINFORMATION.out.bam )

    SAMTOOLS_INDEX( PICARD_ADDORREPLACEREADGROUPS.out.bam )

    FASTQC_POST( FAQCS.out.reads )

    QUALIMAP_BAMQC(
        PICARD_ADDORREPLACEREADGROUPS.out.bam,
        []
    )

    ch_alignment_combined = PICARD_ADDORREPLACEREADGROUPS.out.bam
        .join( SAMTOOLS_INDEX.out.bai )

    SAMTOOLS_STATS    ( ch_alignment_combined, ref_fasta )

    SAMTOOLS_FLAGSTAT ( ch_alignment_combined )

    SAMTOOLS_IDXSTATS ( ch_alignment_combined )

    // Use the full TXT bundle from FAQCS to ensure all required qa.* inputs are staged for QC_REPORT
    ch_qcreport_input = FAQCS.out.txt.join(QUALIMAP_BAMQC.out.results)

    QC_REPORT(ch_qcreport_input, ref_fasta.map{ meta, fasta -> fasta })

    ch_versions = ch_versions.mix(
        SEQKIT_PAIR.out.versions,
        FAQCS.out.versions,
        BWA_MEM.out.versions,
        PICARD_MARKDUPLICATES.out.versions,
        PICARD_CLEANSAM.out.versions,
        PICARD_FIXMATEINFORMATION.out.versions,
        PICARD_ADDORREPLACEREADGROUPS.out.versions,
        SAMTOOLS_INDEX.out.versions,
        FASTQC_POST.out.versions,
        SAMTOOLS_STATS.out.versions,
        QUALIMAP_BAMQC.out.versions,
        QC_REPORT.out.versions,
        SAMTOOLS_IDXSTATS.out.versions,
        SAMTOOLS_FLAGSTAT.out.versions
    )

    if (params.coverage != 0) {
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
    }

    ch_alignment          = PICARD_ADDORREPLACEREADGROUPS.out.bam
    ch_alignment_index    = SAMTOOLS_INDEX.out.bai
    ch_qcreport           = QC_REPORT.out.qc_line

    emit:
    alignment          = ch_alignment                    // channel: [ val(meta), bam ]
    alignment_index    = ch_alignment_index              // channel: [ val(meta), bai ]
    alignment_combined = ch_alignment_combined           // channel: [ val(meta), bam, bai ]
    qualimap           = QUALIMAP_BAMQC.out.results      // channel: [ val(meta), results ]
    stats              = SAMTOOLS_STATS.out.stats        // channel: [ val(meta), stats ]
    flagstat           = SAMTOOLS_FLAGSTAT.out.flagstat  // channel: [ val(meta), flagstat ]
    idxstats           = SAMTOOLS_IDXSTATS.out.idxstats  // channel: [ val(meta), idxstats ]
    post_qc            = FASTQC_POST.out.zip             // channel: [ val(meta), zip ]
    versions           = ch_versions                     // channel: [ ch_versions ]
    qc_lines           = ch_qcreport                     // channel: [ qc_line ]

}
