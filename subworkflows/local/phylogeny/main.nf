//
// Phylogenies subworkflow
//

include { RAPIDNJ          } from '../../../modules/nf-core/rapidnj/main'
include { FASTTREE         } from '../../../modules/nf-core/fasttree/main'
include { IQTREE           } from '../../../modules/nf-core/iqtree/main'
include { RAXMLNG          } from '../../../modules/nf-core/raxmlng/main'
include { QUICKSNP         } from '../../../modules/local/quicksnp/main'
include { MICROREACTSHAPES } from '../../../modules/local/microreactshapes/main'

workflow CREATE_PHYLOGENY {
    take:
    fasta                 // channel: aligned pseudogenomes or filtered version
    constant_sites_string // val: string of constant sites A,C,G,T
    snpdists_tsv          // channel: tsv of SNP distances
    rapidnj               // val: whether to run rapidnj
    fasttree              // val: whether to run fasttree
    iqtree                // val: whether to run iqtree
    raxmlng               // val: whether to run raxmlng
    amdp                  // val: whether to run amdp
    metadata_csv          // val: metadata CSV file for microreact shapes
    geolocation_csv       // val: geolocation CSV file for microreact shapes


    main:
    ch_versions = Channel.empty()

    rapidnj_tree = Channel.empty()
    if (rapidnj) {
        RAPIDNJ (
            fasta
        )
        ch_versions = ch_versions.mix(RAPIDNJ.out.versions)

        rapidnj_tree = RAPIDNJ.out.phylogeny
    }

    fasttree_tree = Channel.empty()
    if (fasttree) {
        FASTTREE (
            fasta
        )
        ch_versions = ch_versions.mix(FASTTREE.out.versions)

        fasttree_tree = FASTTREE.out.phylogeny
    }

    iqtree_tree = Channel.empty()
    if (iqtree) {
        IQTREE (
            [ [], fasta, [] ], [], [], [], [], [], [], [], [], [], [], [], []
        )
        ch_versions = ch_versions.mix(IQTREE.out.versions)

        iqtree_tree = IQTREE.out.phylogeny
    }

    raxmlng_tree = Channel.empty()
    if (raxmlng) {
        RAXMLNG (
            [ [],fasta, 'GTR+G' ]
        )
        ch_versions = ch_versions.mix(RAXMLNG.out.versions)

        raxmlng_tree = RAXMLNG.out.phylogeny
    }

    QUICKSNP (
        snpdists_tsv
    )
    ch_versions = ch_versions.mix(QUICKSNP.out.versions)

    shapes = Channel.empty()
    if (amdp) {
        MICROREACTSHAPES (
            QUICKSNP.out.quicksnp_tree,
            metadata_csv,
            geolocation_csv
        )
        ch_versions = ch_versions.mix(MICROREACTSHAPES.out.versions)

        shapes = MICROREACTSHAPES.out.shapes
    }

    emit:
    rapidnj_tree  = rapidnj_tree                  // channel: [ phylogeny ]
    fasttree_tree = fasttree_tree                 // channel: [ phylogeny ]
    iqtree_tree   = iqtree_tree                   // channel: [ phylogeny ]
    raxmlng_tree  = raxmlng_tree                  // channel: [ phylogeny ]
    quicksnp_tree = QUICKSNP.out.quicksnp_tree    // channel: [ phylogeny ]
    shapes        = shapes                        // channel: [ phylogeny  ]
    versions      = ch_versions                   // channel: [ ch_versions ]
}
