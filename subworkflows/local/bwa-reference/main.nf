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

process INPUT_PROC {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(fasta)

    output:
    tuple val(meta), path("reference.fasta", includeInputs: true), path("reference.copy.fasta"), emit: ref_fasta

    script:
    meta = [
            "id": 'reference'
    ]
    def is_compressed = false
    if(fasta.getName().endsWith(".gz"))
    {
        is_compressed = true
    }
    """
    echo ${meta.id}
    if [[ ${is_compressed} == "true" ]]; then
        gunzip -c $fasta > reference.fasta
    else
        mv ${fasta} reference.fasta
    fi
    cp reference.fasta reference.copy.fasta
    """

}

workflow BWA_REFERENCE {

    take:
    fasta // channel: fasta

    main:

    ch_use_fasta          = Channel.empty()
    ch_masked_fasta       = Channel.empty()
    ch_samtools_index     = Channel.empty()
    ch_bwa_index          = Channel.empty()
    ch_dict               = Channel.empty()
    ch_reference_combined = Channel.empty()
    ch_versions           = Channel.empty()

    INPUT_PROC( fasta )
    NUCMER( INPUT_PROC.out )
    COORDSTOBED( NUCMER.out.delta )
    // BEDTOOLS_MASKFASTA expects tuple val(meta), path(bed) and a path-only fasta
    // Broadcast fasta path alongside each bed entry to ensure pairing
    def ref_fasta_path = INPUT_PROC.out.map{ meta, fa, fa2-> fa }
    def bed_with_fasta = COORDSTOBED.out.bed.map{ meta, bed -> [ meta, bed ] }
    BEDTOOLS_MASKFASTA( bed_with_fasta, ref_fasta_path )


    if(params.mask)
    {
        // else use nucmer masked fasta input
        ch_use_fasta = BEDTOOLS_MASKFASTA.out.fasta
    } else
    {
        // If no_mask is set, use original fasta input
        ch_use_fasta = INPUT_PROC.out.map{meta, fa1, fa2 -> [ meta, fa1 ] }
    }

    BWA_INDEX(ch_use_fasta)
    SAMTOOLS_FAIDX(ch_use_fasta, [[],[]], false)
    PICARD_CREATESEQUENCEDICTIONARY(ch_use_fasta)


    // reference_fasta, samtools_faidx, bwa_index, dict
    // Use join to combine channels by meta key instead of combine
    ch_use_fasta
        .join(SAMTOOLS_FAIDX.out.fai, by: 0)
        .join(BWA_INDEX.out.index, by: 0)
        .join(PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict, by: 0)
        .map { meta, fa1, fai, bai, dict -> [meta, fa1, fai, bai, dict] }
        .set { ch_reference_combined }


    // Collect versions information
    ch_versions           = ch_versions.mix( NUCMER.out.versions,
                                             BWA_INDEX.out.versions,
                                             SAMTOOLS_FAIDX.out.versions,
                                             PICARD_CREATESEQUENCEDICTIONARY.out.versions,
                                             COORDSTOBED.out.versions,
                                             BEDTOOLS_MASKFASTA.out.versions )
    ch_masked_fasta       = ch_use_fasta
    ch_samtools_index     = SAMTOOLS_FAIDX.out.fai
    ch_bwa_index          = BWA_INDEX.out.index
    ch_dict               = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict

    emit:
    masked_fasta       = ch_masked_fasta        // channel: [ val(meta), fas ]
    samtools_index     = ch_samtools_index      // channel: [ val(meta), fai ]
    bwa_index          = ch_bwa_index           // channel: [ val(meta), bwa ]
    dict               = ch_dict                // channel: [ val(meta), dict ]
    reference_combined = ch_reference_combined  // channel: [ val(meta), fa, fai, bai, dict ]
    versions           = ch_versions            // channel: [ ch_versions ]
}
