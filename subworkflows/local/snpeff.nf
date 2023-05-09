//
// Run Snpeff
//

include { SNPEFF as SNPEFF_ANN             } from '../../modules/local/snpeff_local'
include { TABIX_BGZIPTABIX                 } from '../../modules/nf-core/modules/nf-core/tabix/bgziptabix/main'
include { SNPEFFR			   } from '../../modules/local/snpeffr'

workflow SNPEFF {
    take:
    vcf                                             // channel: [val(meta), [ vcf ] ]
    species

    main:

    ch_versions = Channel.empty()

    //
    //SNPEFF
    //

    SNPEFF_ANN (
        vcf,
        species
    )
    ch_snpeff_vcf    = SNPEFF_ANN.out.vcf
    ch_snpeff_csv    = SNPEFF_ANN.out.report
    ch_snpeff_txt    = SNPEFF_ANN.out.genes_txt
    ch_snpeff_html   = SNPEFF_ANN.out.summary_html
    ch_versions      = ch_versions.mix(SNPEFF_ANN.out.versions.first())

    //
    //ZIP & INDEXING VCF
    //

    TABIX_BGZIPTABIX (
        ch_snpeff_vcf
    )
    ch_tabix_tbi     = TABIX_BGZIPTABIX.out.gz_tbi
    ch_tabix_vcf     = TABIX_BGZIPTABIX.out.gz
    ch_versions      = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    SNPEFFR (
        ch_tabix_vcf
    )
    
    ch_snpeffr_csv    = SNPEFFR.out.report
    ch_versions       = ch_versions.mix(SNPEFFR.out.versions.first())

    emit:
    csv              = ch_snpeff_csv                   // channel: [ val(meta), [ csv ] ]
    csv_snpeffr      = ch_snpeffr_csv                  // channel: [ val(meta), [ csv ] ]
    txt              = ch_snpeff_txt                   // channel: [ val(meta), [ txt ] ]
    html             = ch_snpeff_html                  // channel: [ val(meta), [ html ] ]
    vcf_tbi          = ch_tabix_tbi                    // channel: [ val(meta), [ tbi ] ]           
    vcf              = ch_tabix_vcf                    // channel: [ val(meta), [ gz ] ]
    versions         = ch_versions                     // channel: [ versions.yml ]

}


