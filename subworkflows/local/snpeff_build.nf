//
// Run Snpeff Build
//

include { SNPEFF_BUILD        } from '../../modules/local/snpeff_build.nf'
include { SNPEFF              } from '../../modules/nf-core/modules/nf-core/snpeff/main'
include { TABIX_BGZIPTABIX    } from '../../modules/nf-core/modules/nf-core/tabix/bgziptabix/main'


workflow SNPEFF_DB {
    take:
    fasta                                           // path   : genome.fasta
    gff                                             // path   : genome.gff

    main:

    ch_versions = Channel.empty()

    //
    //SNPEFF
    //
    
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    SNPEFF_BUILD (
        gff,
        fasta
    )
    ch_snpeff_db     = SNPEFF_BUILD.out.bin
    ch_snpeff_config = SNPEFF_BUILD.out.config
    ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)

    //
    //SNPEFF
    //
    ch_snpeff_vcf    = Channel.empty()
    ch_snpeff_csv    = Channel.empty()
    ch_snpeff_txt    = Channel.empty()
    ch_snpeff_html   = Channel.empty()
    SNPEFF (
        ch_snpeff_db,
        ch_snpeff_config,
        vcf,
        fasta // needed ?
    )
    ch_snpeff_vcf    = SNPEFF.out.vcf
    ch_snpeff_csv    = SNPEFF.out.csv
    ch_snpeff_txt    = SNPEFF.out.txt
    ch_snpeff_html   = SNPEFF.out.html
    ch_versions      = ch_versions.mix(SNPEFF.out.versions.first())

    //
    //ZIP & INDEXING VCF
    //
    ch_tabix_tbi     = Channel.empty()
    ch_tabix_vcf     = Channel.empty()
    TABIX_BGZIPTABIX (
        ch_snpeff_vcf
    )
    ch_tabix_tbi     = TABIX_BGZIPTABIX.out.gz_tbi
    ch_tabix_vcf     = TABIX_BGZIPTABIX.out.gz
    ch_versions      = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())


    emit:         
    snpeff_db        = ch_snpeff_db                    // channel: [ val(meta), [ csv ] ]
    config           = ch_snpeff_config                // channel: [ val(meta), [ txt ] ]
    csv              = ch_snpeff_csv                   // channel: [ val(meta), [ csv ] ]
    txt              = ch_snpeff_txt                   // channel: [ val(meta), [ txt ] ]
    html             = ch_snpeff_html                  // channel: [ val(meta), [ html ] ]
    vcf_tbi          = ch_tabix_tbi                    // channel: [ val(meta), [ tbi ] ]           
    vcf              = ch_tabix_vcf                    // channel: [ val(meta), [ vcf.gz ] ]
    versions         = ch_versions                     // channel: [ versions.yml ]

}

