/*
========================================================================================
    BWA Reference Sub-Workflow
========================================================================================
*/

include { NUCMER }                          from '../../modules/nf-core/modules/nucmer/main'
// include { SHOW_COORDS }                     from '../../modules/local/showcoords.nf'
// include { COORDSTOBED }                     from '../../modules/local/coordstobed.nf'
// include { BEDTOOLS_MASKFASTA }              from '../../modules/nf-core/modules/bedtools/maskfasta/'
include { BWA_INDEX }                       from '../../modules/nf-core/modules/bwa/index/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/modules/picard/createsequencedictionary/main'
include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/modules/samtools/faidx/main'

process INPUT_PROC {

    input:
    path(fasta)

    output:
    tuple val(meta), path(fasta, includeInputs: true), path("*.copy.fasta")

    script:
    meta = [
            "id": 'reference'
    ]

    """
    echo ${meta.id}
    cp ${fasta} ${fasta}.copy.fasta
    """

}

workflow BWA_REFERENCE {

    take:
    fasta

    main:
    ch_versions = Channel.empty()

    INPUT_PROC(fasta)
    NUCMER( INPUT_PROC.out )
    
    // TODO: Add masking and use that for BWA_INDEX, SAMTOOLS, PICARD
    // SHOW_COORDS()
    // COORDSTOBED()
    // BEDTOOLS_MASKFASTA()

    BWA_INDEX(fasta)
    SAMTOOLS_FAIDX(INPUT_PROC.out.map{meta, fa1, fa2->[meta, fa1] })
    PICARD_CREATESEQUENCEDICTIONARY(INPUT_PROC.out.map{meta, fa1, fa2->[meta, fa1] })

    // Collect versions information
    ch_versions = ch_versions.mix(  NUCMER.out.versions, 
                                    BWA_INDEX.out.versions, 
                                    SAMTOOLS_FAIDX.out.versions, 
                                    PICARD_CREATESEQUENCEDICTIONARY.out.versions
                                )
    

    emit:
    samtools_index = SAMTOOLS_FAIDX.out.fai
    bwa_index = BWA_INDEX.out.index
    dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
    versions = ch_versions // channel: [ versions.yml ]
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


