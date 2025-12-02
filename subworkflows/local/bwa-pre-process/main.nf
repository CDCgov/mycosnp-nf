/*
========================================================================================
    BWA Pre-Process Sub-Workflow
========================================================================================
*/

include { SEQKIT_PAIR                          } from '../../../modules/nf-core/seqkit/pair/main'
include { SEQTK_SAMPLE                         } from '../../../modules/nf-core/seqtk/sample/main'
include { FAQCS                                } from '../../../modules/nf-core/faqcs/main'
include { BWA_MEM                              } from '../../../modules/nf-core/bwa/mem/main'
include { PICARD_MARKDUPLICATES                } from '../../../modules/nf-core/picard/markduplicates/main'
include { PICARD_CLEANSAM                      } from '../../../modules/nf-core/picard/cleansam/main'
include { PICARD_FIXMATEINFORMATION            } from '../../../modules/nf-core/picard/fixmateinformation/main'
include { PICARD_ADDORREPLACEREADGROUPS        } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { FASTQC as FASTQC_POST                } from '../../../modules/nf-core/fastqc/main'
include { DOWNSAMPLE_RATE                      } from '../../../modules/local/downsample_rate/main'
include { SAMTOOLS_STATS                       } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS                    } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT                    } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_COVERAGE                    } from '../../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_DEPTH                       } from '../../../modules/nf-core/samtools/depth/main'
include { QC_REPORT                            } from '../../../modules/local/qc_report/main'


// Function to pair each item in a channel with a reference file
def pairWithReference(input_channel, reference_channel) {
    return input_channel
        .map { meta, _file -> meta }
        .combine( reference_channel.map { _m, ref -> ref } )
        .map { meta, ref -> tuple(meta, ref) }
}

/*
========================================================================================
    Workflow
========================================================================================
*/

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
    DOWNSAMPLE_RATE(
        SEQKIT_PAIR.out.reads,
        ref_fasta.map{ _meta, fa -> fa },
        params.coverage
    )

    ch_versions = ch_versions.mix(DOWNSAMPLE_RATE.out.versions)

    ch_seq_samplerate = SEQKIT_PAIR.out.reads.join(
        DOWNSAMPLE_RATE.out.downsampled_rate.map{
            meta, _sr, snr -> [ meta, snr]
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
    // Reference index and fasta for nf-core BWA_MEM (requires reads,index,fasta,sort flag)
    ch_bwa_index = pairWithReference(FAQCS.out.reads, ref_bwa )

    ch_bwa_fasta = pairWithReference(FAQCS.out.reads, ref_fasta )

    BWA_MEM(
        FAQCS.out.reads,
        ch_bwa_index,
        ch_bwa_fasta,
        true
    )

    ch_versions = ch_versions.mix(
        BWA_MEM.out.versions
    )

    // Reference fasta / fai for MarkDuplicates (nf-core variant requires both)
    ch_md_ref_fasta = pairWithReference(BWA_MEM.out.bam, ref_fasta)
    ch_md_ref_fai = pairWithReference(BWA_MEM.out.bam, ref_fai)

    PICARD_MARKDUPLICATES(
        BWA_MEM.out.bam,
        ch_md_ref_fasta,
        ch_md_ref_fai
    )

    ch_versions = ch_versions.mix( PICARD_MARKDUPLICATES.out.versions )

    PICARD_CLEANSAM (
        PICARD_MARKDUPLICATES.out.bam
    )

    ch_versions = ch_versions.mix(PICARD_CLEANSAM.out.versions)

    PICARD_FIXMATEINFORMATION (
        PICARD_CLEANSAM.out.bam
    )

    ch_versions = ch_versions.mix(PICARD_FIXMATEINFORMATION.out.versions)

    // Reference fasta / fai so each BAM gets paired with singleton reference files
    ch_ref_fasta = pairWithReference(PICARD_FIXMATEINFORMATION.out.bam, ref_fasta)
    ch_ref_fai = pairWithReference(PICARD_FIXMATEINFORMATION.out.bam, ref_fai)

    // nf-core module expects: bam, fasta, fasta_index
    PICARD_ADDORREPLACEREADGROUPS(
        PICARD_FIXMATEINFORMATION.out.bam,
        ch_ref_fasta,
        ch_ref_fai
    )

    ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions)

    FASTQC_POST (
        FAQCS.out.reads
    )

    ch_versions = ch_versions.mix(FASTQC_POST.out.versions)

    ch_alignment_combined = PICARD_ADDORREPLACEREADGROUPS.out.bam
        .join( PICARD_ADDORREPLACEREADGROUPS.out.bai )

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

    // SAMTOOLS_COVERAGE needs bam, bai, fasta, and fai
    ch_coverage_fasta = pairWithReference(PICARD_ADDORREPLACEREADGROUPS.out.bam, ref_fasta)
    ch_coverage_fai = pairWithReference(PICARD_ADDORREPLACEREADGROUPS.out.bam, ref_fai)

    // Create combined channel with bam+bai for SAMTOOLS_COVERAGE
    ch_coverage_combined = PICARD_ADDORREPLACEREADGROUPS.out.bam
        .join( PICARD_ADDORREPLACEREADGROUPS.out.bai )

    SAMTOOLS_COVERAGE (
        ch_coverage_combined,
        ch_coverage_fasta,
        ch_coverage_fai
    )

    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)

    // SAMTOOLS_DEPTH needs bam and optional intervals
    // Pass bam files and a matching channel with null intervals for each sample
    ch_depth_intervals = PICARD_ADDORREPLACEREADGROUPS.out.bam
        .map { meta, bam -> tuple(meta, []) }

    SAMTOOLS_DEPTH (
        PICARD_ADDORREPLACEREADGROUPS.out.bam,
        ch_depth_intervals
    )

    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

    // QC_REPORT needs FAQCS stats, debug, samtools coverage, samtools depth, and samtools stats
    ch_qcreport_input = FAQCS.out.stats
        .join(FAQCS.out.debug)
        .join(SAMTOOLS_COVERAGE.out.coverage)
        .join(SAMTOOLS_DEPTH.out.tsv)
        .join(SAMTOOLS_STATS.out.stats)

    QC_REPORT(
        ch_qcreport_input,
        ref_fasta.map{ _meta, fasta -> fasta }
    )

    ch_versions = ch_versions.mix(
        SEQKIT_PAIR.out.versions,
        FAQCS.out.versions,
        BWA_MEM.out.versions,
        PICARD_MARKDUPLICATES.out.versions,
        PICARD_CLEANSAM.out.versions,
        PICARD_FIXMATEINFORMATION.out.versions,
        PICARD_ADDORREPLACEREADGROUPS.out.versions,
        FASTQC_POST.out.versions,
        SAMTOOLS_STATS.out.versions,
        QC_REPORT.out.versions,
        SAMTOOLS_IDXSTATS.out.versions,
        SAMTOOLS_FLAGSTAT.out.versions,
        SAMTOOLS_COVERAGE.out.versions,
        SAMTOOLS_DEPTH.out.versions
    )

    if (params.coverage != 0) {
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
    }

    ch_alignment          = PICARD_ADDORREPLACEREADGROUPS.out.bam
    ch_alignment_index    = PICARD_ADDORREPLACEREADGROUPS.out.bai
    ch_qcreport           = QC_REPORT.out.qc_line

    emit:
    alignment          = PICARD_ADDORREPLACEREADGROUPS.out.bam  // channel: [ val(meta), bam ]
    alignment_combined = ch_alignment_combined                  // channel: [ val(meta), bam, bai ]
    stats              = SAMTOOLS_STATS.out.stats               // channel: [ val(meta), stats ]
    flagstat           = SAMTOOLS_FLAGSTAT.out.flagstat         // channel: [ val(meta), flagstat ]
    idxstats           = SAMTOOLS_IDXSTATS.out.idxstats         // channel: [ val(meta), idxstats ]
    coverage           = SAMTOOLS_COVERAGE.out.coverage         // channel: [ val(meta), txt ]
    depth              = SAMTOOLS_DEPTH.out.tsv                 // channel: [ val(meta), tsv ]
    post_qc            = FASTQC_POST.out.zip                    // channel: [ val(meta), zip ]
    versions           = ch_versions                            // channel: [ ch_versions ]
    qc_lines           = QC_REPORT.out.qc_line                  // channel: [ qc_line ]

}
