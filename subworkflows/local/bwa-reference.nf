/*
========================================================================================
    BWA Reference Sub-Workflow
========================================================================================
*/

include { NUCMER }                          from '../../modules/nf-core/modules/nucmer/main'
include { SHOW_COORDS }                     from '../../modules/local/showcoords.nf'
// include { COORDSTOBED }                     from '../../modules/local/coordstobed.nf'
// include { BEDTOOLS_MASKFASTA }              from '../../modules/nf-core/modules/bedtools/maskfasta/'
include { BWA_INDEX }                       from '../../modules/nf-core/modules/bwa/index/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/modules/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/modules/samtools/faidx/main'

process INPUT_PROC {

    input:
    path(fasta)

    output:
    tuple val(meta), path("reference.fasta", includeInputs: true), path("reference.copy.fasta")

    script:
    meta = [
            "id": 'reference'
    ]

    """
    echo ${meta.id}
    mv ${fasta} reference.fasta
    cp reference.fasta reference.copy.fasta
    """

}

workflow BWA_REFERENCE {

    take:
    fasta // channel: [ val(meta), [ fastq ] ]

    main:
    
    ch_masked_fasta       = Channel.empty()
    ch_samtools_index     = Channel.empty()
    ch_bwa_index          = Channel.empty()
    ch_dict               = Channel.empty()
    ch_reference_combined = Channel.empty()
    ch_versions           = Channel.empty()

    INPUT_PROC(fasta)
    NUCMER( INPUT_PROC.out )
    // NUCMER.out.coords
    // NUCMER.out.delta

    // TODO: Add masking and use that for BWA_INDEX, SAMTOOLS, PICARD
    SHOW_COORDS(NUCMER.out.delta)
    // COORDSTOBED()
    // BEDTOOLS_MASKFASTA()

    BWA_INDEX(INPUT_PROC.out.map{meta, fa1, fa2 -> fa1})
    SAMTOOLS_FAIDX(INPUT_PROC.out.map{meta, fa1, fa2->[meta, fa1] })
    PICARD_CREATESEQUENCEDICTIONARY(INPUT_PROC.out.map{meta, fa1, fa2->[meta, fa1] })

    
    // reference_fasta, samtools_faidx, bwa_index, dict
    INPUT_PROC.out.combine(SAMTOOLS_FAIDX.out.fai).combine(BWA_INDEX.out.index).combine(PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict)
      .map{meta, fa1, fa2, meta2, fai, bai, meta4, dict -> [meta, fa1, fai, bai, dict] }
      .set{ch_reference_combined}
      
    // ch_combined.view()


    // Collect versions information
    ch_versions           = ch_versions.mix( NUCMER.out.versions, 
                                             BWA_INDEX.out.versions, 
                                             SAMTOOLS_FAIDX.out.versions, 
                                             PICARD_CREATESEQUENCEDICTIONARY.out.versions )
    ch_masked_fasta       = fasta
    ch_samtools_index     = SAMTOOLS_FAIDX.out.fai 
    ch_bwa_index          = BWA_INDEX.out.index 
    ch_dict               = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict

    emit:
    masked_fasta       = ch_masked_fasta        // channel: [ val(meta), [ fai ] ]
    samtools_index     = ch_samtools_index      // channel: [ val(meta), [ fai ] ]
    bwa_index          = ch_bwa_index           // channel: [ val(meta), [ index ] ]
    dict               = ch_dict                // channel: [ val(meta), [ dict ] ]
    reference_combined = ch_reference_combined  // channel: [ val(meta), [ vcf ] ]
    versions           = ch_versions            // channel: [ ch_versions ]
}


/*
## Step 1: Mask Repeats
  mask-repeats:
    git: https://gitlab.com/geneflow/apps/mask-repeats-gf2.git
    version: '0.2'
nucmer -p ${OUTPUT_BASE} --maxmatch --nosimplify -t $THREADS $REFERENCE_SEQUENCE $REFERENCE_SEQUENCE

## Step 2: Index the masked Reference Sequence with Picard

  index-reference:
    git: https://gitlab.com/geneflow/apps/index-reference-gf2.git
    version: '0.1'
picard CreateSequenceDictionary R=${OUTPUT_BASE}.fasta O=${OUTPUT_BASE}.dict

---NOTE: THIS SHOULD BE AN OPTIONAL STEP - NEED OPTION TO SKIP MASKING
show-coords -r -T $OUTPUT_BASE.delta -H
awk '{if (\$1 != \$3 && \$2 != \$4) print \$0}'
awk '{print \$8\"\\t\"\$1\"\\t\"\$2}'
STDOUT  ->  ${OUTPUT_FULL}/${OUTPUT_BASE}.bed
bedtools maskfasta -fi $REFERENCE_SEQUENCE -bed $BEDFILE.bed -fo $NEWOUTPUT.fasta
--- END OPTIONAL STEP

## Step 3: Index masked Reference sequence with samtools faidx
samtools faidx ${OUTPUT_BASE}.fasta

## Step 4: Create BWA index
  bwa-index:
    git: https://gitlab.com/geneflow/apps/bwa-index-gf2.git
    version: '0.7.17-03'
bwa index -p $OUTPUT_BASE $REFERENCE

*/


