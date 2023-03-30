//
// Phylogenies subworkflow
//

include { RAPIDNJ  } from '../../modules/nf-core/modules/rapidnj/main'
include { FASTTREE } from '../../modules/nf-core/modules/fasttree/main'
include { IQTREE   } from '../../modules/nf-core/modules/iqtree/main'
include { RAXMLNG  } from '../../modules/nf-core/modules/raxmlng/main'
include { QUICKSNP } from '../../modules/local/quicksnp.nf'

workflow CREATE_PHYLOGENY {
    take:
    fasta                 // channel: aligned pseudogenomes or filtered version
    constant_sites_string // val: string of constant sites A,C,G,T
    snpdists_tsv

    main:
    ch_versions = Channel.empty()

    rapidnj_tree = Channel.empty()

    if (params.rapidnj) {
        RAPIDNJ(fasta)
        rapidnj_tree = RAPIDNJ.out.phylogeny
        ch_versions = ch_versions.mix(RAPIDNJ.out.versions)
    }

    fasttree_tree = Channel.empty()

    if (params.fasttree) {
        FASTTREE(fasta)
        fasttree_tree = FASTTREE.out.phylogeny
        ch_versions = ch_versions.mix(FASTTREE.out.versions)
    }

    iqtree_tree = Channel.empty()
    if (params.iqtree) {
        IQTREE(fasta, constant_sites_string)
        iqtree_tree = IQTREE.out.phylogeny
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    raxmlng_tree = Channel.empty()
    if (params.raxmlng) {
        RAXMLNG(fasta)
        raxmlng_tree = RAXMLNG.out.phylogeny
        ch_versions = ch_versions.mix(RAXMLNG.out.versions)
    }

    QUICKSNP(snpdists_tsv)

    emit:
    rapidnj_tree      = rapidnj_tree     // channel: [ phylogeny ]
    fasttree_tree     = fasttree_tree    // channel: [ phylogeny ]
    iqtree_tree       = iqtree_tree      // channel: [ phylogeny ]
    raxmlng_tree      = raxmlng_tree     // channel: [ phylogeny ]
    //quicksnp_tree     = quicksnp_tree    // channel: [ phylogeny ]
    versions          = ch_versions      // channel: [ ch_versions ]
}
