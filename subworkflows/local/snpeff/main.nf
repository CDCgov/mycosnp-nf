//
// Run Snpeff
//

include { SNPEFF_SNPEFF as SNPEFF_ANN      } from '../../../modules/local/snpeff/snpeff/main'
include { TABIX_BGZIPTABIX                 } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { SNPEFFR			               } from '../../../modules/local/snpeffr/main'

workflow SNPEFF {
    take:
    vcf                                             // channel: [val(meta), [ vcf ] ]
    species                                         // val: species name for SnpEff
    snpeffcache                                     // channel: path to snpeff cache directory
    ref_name                                        // val: reference name for output naming in snpeffr

    main:
    ch_versions = Channel.empty()

    SNPEFF_ANN (
        vcf,
        species,
        snpeffcache.map{ cache -> [[:], cache] }.first()
    )
    ch_versions = ch_versions.mix(SNPEFF_ANN.out.versions)

    //ZIP & INDEXING VCF
    TABIX_BGZIPTABIX (
        SNPEFF_ANN.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    ch_tabix_vcf = TABIX_BGZIPTABIX.out.gz_tbi
        .mix( TABIX_BGZIPTABIX.out.gz_csi )
        .map { meta, gz, idx -> [ meta, gz ] }

    SNPEFFR (
        ch_tabix_vcf,
        ref_name
    )
    ch_versions = ch_versions.mix(SNPEFFR.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
