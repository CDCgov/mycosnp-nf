/*
========================================================================================
    BWA Pre-Process Sub-Workflow
========================================================================================
*/

include { SEQKIT_PAIR                          } from '../../modules/nf-core/modules/seqkit/pair/main'
include { SEQTK_SAMPLE                         } from '../../modules/local/seqtk_sample.nf'
include { FAQCS                                } from '../../modules/nf-core/modules/faqcs/main'
include { BWA_INDEX                            } from '../../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM                              } from '../../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_SORT                        } from '../../modules/nf-core/modules/samtools/sort/main'
include { PICARD_MARKDUPLICATES                } from '../../modules/nf-core/modules/picard/markduplicates/main'
include { PICARD_CLEANSAM                      } from '../../modules/nf-core/modules/picard/cleansam/main'
include { SAMTOOLS_VIEW as PICARDDUPTOCLEANSAM } from '../../modules/nf-core/modules/samtools/view/main'
include { PICARD_FIXMATEINFORMATION            } from '../../modules/nf-core/modules/picard/fixmateinformation/main'
include { PICARD_ADDORREPLACEREADGROUPS        } from '../../modules/nf-core/modules/picard/addorreplacereadgroups/main'
include { SAMTOOLS_INDEX                       } from '../../modules/nf-core/modules/samtools/index/main'
include { FASTQC as FASTQC_POST                } from '../../modules/nf-core/modules/fastqc/main'
include { QUALIMAP_BAMQC                       } from '../../modules/nf-core/modules/qualimap/bamqc/main'
include { DOWNSAMPLE_RATE                      } from '../../modules/local/downsample_rate.nf'
include { SAMTOOLS_STATS                       } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS                    } from '../../modules/nf-core/modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT                    } from '../../modules/nf-core/modules/samtools/flagstat/main'
include { QC_REPORT                            } from '../../modules/local/qc_report.nf'



workflow BWA_PREPROCESS {

    take:
    reference //channel: [ tuple reference_fasta, samtools_faidx, bwa_index ]
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
    DOWNSAMPLE_RATE(SEQKIT_PAIR.out.reads, reference[0], params.coverage)
    ch_seq_samplerate = SEQKIT_PAIR.out.reads.join(
        DOWNSAMPLE_RATE.out.downsampled_rate.map{ meta, sr, snr -> [ meta, snr]}
        )
    SEQTK_SAMPLE(ch_seq_samplerate)
    FAQCS(SEQTK_SAMPLE.out.reads)
    }
    
    BWA_MEM(FAQCS.out.reads, reference[2], true)
    PICARD_MARKDUPLICATES(BWA_MEM.out.bam)
    PICARD_CLEANSAM(PICARD_MARKDUPLICATES.out.bam)
    PICARD_FIXMATEINFORMATION(PICARD_CLEANSAM.out.bam)
    PICARD_ADDORREPLACEREADGROUPS(PICARD_FIXMATEINFORMATION.out.bam)
    SAMTOOLS_INDEX(PICARD_ADDORREPLACEREADGROUPS.out.bam)
    FASTQC_POST(FAQCS.out.reads)
    QUALIMAP_BAMQC(PICARD_ADDORREPLACEREADGROUPS.out.bam, [], false)

    ch_alignment_combined = PICARD_ADDORREPLACEREADGROUPS.out.bam.join(SAMTOOLS_INDEX.out.bai)
    
    SAMTOOLS_STATS    ( ch_alignment_combined, reference[0] )
    SAMTOOLS_FLAGSTAT ( ch_alignment_combined )
    SAMTOOLS_IDXSTATS ( ch_alignment_combined )
    
    ch_qcreport_input = FAQCS.out.txt.join(QUALIMAP_BAMQC.out.results)

    QC_REPORT(ch_qcreport_input, reference[0])

    ch_versions            = ch_versions.mix(  SEQKIT_PAIR.out.versions, 
                                               FAQCS.out.versions,
                                               BWA_MEM.out.versions,
                                               PICARD_MARKDUPLICATES.out.versions,
                                               PICARD_CLEANSAM.out.versions,
                                               PICARD_FIXMATEINFORMATION.out.versions,
                                               PICARD_ADDORREPLACEREADGROUPS.out.versions,
                                               SAMTOOLS_INDEX.out.versions,
                                               FASTQC_POST.out.versions,
                                               SAMTOOLS_STATS.out.versions,
                                               QUALIMAP_BAMQC.out.versions
                                            )
    if (params.coverage != 0) { 
        ch_versions.mix(SEQTK_SAMPLE.out.versions)
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
