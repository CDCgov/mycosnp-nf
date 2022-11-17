//
// Run Snpeff
//

include { SNPEFF              } from '../../modules/nf-core/modules/nf-core/snpeff/main'
include { TABIX_BGZIPTABIX    } from '../../modules/nf-core/modules/nf-core/tabix/bgziptabix/main'

workflow SNPEFF {
    take:
    vcf                                             // channel: [val(meta), [ vcf ] ]
    snpeffdb                                        // path   : snpEff database
    snpeffconfig                                    // path   : snpEff config

    main:

    ch_versions = Channel.empty()

    //
    //SNPEFF
    //
    println(vcf)
    println(snpeffdb)
    println(snpeffconfig)

    
    SNPEFF (
        vcf,
        snpeffdb,
        snpeffconfig
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
    csv              = ch_snpeff_csv                   // channel: [ val(meta), [ csv ] ]
    txt              = ch_snpeff_txt                   // channel: [ val(meta), [ txt ] ]
    html             = ch_snpeff_html                  // channel: [ val(meta), [ html ] ]
    vcf_tbi          = ch_tabix_tbi                    // channel: [ val(meta), [ tbi ] ]           
    vcf              = ch_tabix_vcf                    // channel: [ val(meta), [ vcf.gz ] ]
    versions         = ch_versions                     // channel: [ versions.yml ]

}

