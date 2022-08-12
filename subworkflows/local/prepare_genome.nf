//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GFF          } from '../../modules/nf-core/modules/gunzip/main'
include { CUSTOM_GETCHROMSIZES          } from '../../modules/nf-core/modules/custom/getchromsizes/main'
include { SNPEFF_BUILD                  } from '../../modules/local/snpeff_build'

workflow PREPARE_GENOME {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            [ [:], params.fasta ]
        )
        ch_fasta    = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GFF annotation file
    //
    if (params.ref_gff) {
        if (params.ref_gff.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], params.ref_gff ]
            )
            ch_gff      = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = file(params.ref_gff)
        }
    } else {
        ch_gff = []
    }

    //
    // Create chromosome sizes file
    //
    ch_fai         = Channel.empty()
    ch_chrom_sizes = Channel.empty()
    CUSTOM_GETCHROMSIZES (
        ch_fasta
    )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Make snpEff database
    //
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (params.ref_gff && !params.skip_snpeff) {
        SNPEFF_BUILD (
            ch_fasta,
            ch_gff
        )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
        ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)
    }

    emit:
    fasta                = ch_fasta                // path: genome.fasta
    gff                  = ch_gff                  // path: genome.gff
    chrom_sizes          = ch_chrom_sizes          // path: genome.sizes
    snpeff_db            = ch_snpeff_db            // path: snpeff_db
    snpeff_config        = ch_snpeff_config        // path: snpeff.config

    versions             = ch_versions             // channel: [ versions.yml ]
}
