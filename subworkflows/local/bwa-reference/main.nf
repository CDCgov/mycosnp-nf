/*
========================================================================================
    BWA Reference Sub-Workflow
========================================================================================
*/

include { NUCMER                          } from '../../../modules/nf-core/nucmer/main'
include { COORDSTOBED                     } from '../../../modules/local/coordstobed/main'
include { BEDTOOLS_MASKFASTA              } from '../../../modules/nf-core/bedtools/maskfasta/main'
include { BWA_INDEX                       } from '../../../modules/nf-core/bwa/index/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                  } from '../../../modules/nf-core/samtools/faidx/main'
include { INPUT_PROC                      } from '../../../modules/local/input_proc/main'

workflow BWA_REFERENCE {

    take:
    fasta // channel: fasta
    mask  // val: whether to mask reference based on nucmer results

    main:
    ch_versions = Channel.empty()

    INPUT_PROC (
        fasta
    )
    ch_versions = ch_versions.mix(INPUT_PROC.out.versions)

    NUCMER (
        INPUT_PROC.out.ref_fasta
    )
    ch_versions = ch_versions.mix(NUCMER.out.versions)

    COORDSTOBED (
        NUCMER.out.delta
    )
    ch_versions = ch_versions.mix(COORDSTOBED.out.versions)

    // Broadcast fasta path alongside each bed entry to ensure pairing
    def ref_fasta_path = INPUT_PROC.out.ref_fasta.map {
        meta, fa, fa2-> fa
    }

    def bed_with_fasta = COORDSTOBED.out.bed.map {
        meta, bed -> [ meta, bed ]
    }

    BEDTOOLS_MASKFASTA (
        bed_with_fasta, ref_fasta_path
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions)

    if(mask) {
        // else use nucmer masked fasta input
        ch_use_fasta = BEDTOOLS_MASKFASTA.out.fasta
    } else {
        // If no_mask is set, use original fasta input
        ch_use_fasta = INPUT_PROC.out.ref_fasta.map {
            meta, fa1, fa2 -> [ meta, fa1 ]
        }
    }

    BWA_INDEX (
        ch_use_fasta
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    SAMTOOLS_FAIDX (
        ch_use_fasta, [[],[]], false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    PICARD_CREATESEQUENCEDICTIONARY (
        ch_use_fasta
    )
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    // Use join to combine channels by meta key instead of combine
    ch_use_fasta
        .join(SAMTOOLS_FAIDX.out.fai)
        .join(BWA_INDEX.out.index)
        .join(PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict)
        .map { meta, fa1, fai, bai, dict -> [meta, fa1, fai, bai, dict] }
        .set { ch_reference_combined }

    emit:
    masked_fasta       = ch_use_fasta                                       // channel: [ val(meta), fas ]
    samtools_index     = SAMTOOLS_FAIDX.out.fai                             // channel: [ val(meta), fai ]
    bwa_index          = BWA_INDEX.out.index                                // channel: [ val(meta), bwa ]
    dict               = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict // channel: [ val(meta), dict ]
    reference_combined = ch_reference_combined                              // channel: [ val(meta), fa, fai, bai, dict ]
    versions           = ch_versions                                        // channel: [ ch_versions ]
}
