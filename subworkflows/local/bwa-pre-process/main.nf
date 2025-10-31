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
    coverage  // val: desired coverage for downsampling (0 = no downsampling)

    main:
    ch_versions           = Channel.empty()

    SEQKIT_PAIR (
        reads
    )
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions)

    if (coverage == 0) {
        FAQCS (
            SEQKIT_PAIR.out.reads
        )
        ch_versions = ch_versions.mix(FAQCS.out.versions)
    }
    else {
        // DOWNSAMPLE_RATE expects a plain path for reference_fasta
        DOWNSAMPLE_RATE (
            SEQKIT_PAIR.out.reads,
            ref_fasta.map{ meta, fa -> fa },
            coverage
        )
        ch_versions = ch_versions.mix(DOWNSAMPLE_RATE.out.versions)

        ch_seq_samplerate = SEQKIT_PAIR.out.reads.join(
            DOWNSAMPLE_RATE.out.downsampled_rate.map {
                meta, sr, snr -> [ meta, snr ]
            }
        )

        SEQTK_SAMPLE(
            ch_seq_samplerate
        )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)

        FAQCS (
            SEQTK_SAMPLE.out.reads
        )
        ch_versions = ch_versions.mix(FAQCS.out.versions)

    }
    // Group trimmed reads by sample so BWA_MEM receives [meta, [R1,R2]]
    // Pass only reads + index path + sort flag as required by local BWA_MEM (3 inputs)
    BWA_MEM (
        FAQCS.out.reads,
        ref_bwa.map{ meta, idx -> idx },
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    PICARD_MARKDUPLICATES (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    PICARD_CLEANSAM (
        PICARD_MARKDUPLICATES.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_CLEANSAM.out.versions)

    PICARD_FIXMATEINFORMATION (
        PICARD_CLEANSAM.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_FIXMATEINFORMATION.out.versions)

    PICARD_ADDORREPLACEREADGROUPS (
        PICARD_FIXMATEINFORMATION.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions)

    SAMTOOLS_INDEX (
        PICARD_ADDORREPLACEREADGROUPS.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    FASTQC_POST (
        FAQCS.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC_POST.out.versions)

    QUALIMAP_BAMQC (
        PICARD_ADDORREPLACEREADGROUPS.out.bam,
        []
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)

    ch_alignment_combined = PICARD_ADDORREPLACEREADGROUPS.out.bam
        .join( SAMTOOLS_INDEX.out.bai )

    SAMTOOLS_STATS (
        ch_alignment_combined,
        ref_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    SAMTOOLS_FLAGSTAT (
        ch_alignment_combined
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    SAMTOOLS_IDXSTATS (
        ch_alignment_combined
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    // Use the full TXT bundle from FAQCS to ensure all required qa.* inputs are staged for QC_REPORT
    ch_qcreport_input = FAQCS.out.txt.join(QUALIMAP_BAMQC.out.results)

    QC_REPORT (
        ch_qcreport_input,
        ref_fasta.map{ meta, fasta -> fasta }
    )
    ch_versions = ch_versions.mix(QC_REPORT.out.versions)

    emit:
    alignment          = PICARD_ADDORREPLACEREADGROUPS.out.bam  // channel: [ val(meta), bam ]
    alignment_index    = SAMTOOLS_INDEX.out.bai                 // channel: [ val(meta), bai ]
    alignment_combined = ch_alignment_combined                  // channel: [ val(meta), bam, bai ]
    qualimap           = QUALIMAP_BAMQC.out.results             // channel: [ val(meta), results ]
    stats              = SAMTOOLS_STATS.out.stats               // channel: [ val(meta), stats ]
    flagstat           = SAMTOOLS_FLAGSTAT.out.flagstat         // channel: [ val(meta), flagstat ]
    idxstats           = SAMTOOLS_IDXSTATS.out.idxstats         // channel: [ val(meta), idxstats ]
    post_qc            = FASTQC_POST.out.zip                    // channel: [ val(meta), zip ]
    versions           = ch_versions                            // channel: [ ch_versions ]
    qc_lines           = QC_REPORT.out.qc_line                  // channel: [ qc_line ]

}
