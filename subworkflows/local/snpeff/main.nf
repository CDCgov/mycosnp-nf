//
// Run Snpeff
//

include { SNPEFF_SNPEFF as SNPEFF_ANN      } from '../../../modules/nf-core/snpeff/snpeff/main'
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

    // Filter VCFs to ensure they contain at least one variant record before running SNPEFFR
    // SNPEFFR fails on empty VCFs (header-only files with no variant data)
    ch_has_variants = ch_tabix_vcf
        .splitText()
        .filter { meta, line -> !line.startsWith('#') }
        .map { meta, line -> [ meta.id, true ] }
        .unique { it[0] }

    ch_tabix_vcf_with_variants = ch_tabix_vcf
        .map { meta, gz -> [ meta.id, meta, gz ] }
        .join( ch_has_variants, by: 0, remainder: true )
        .filter { id, meta, gz, hasVariants ->
            if( !hasVariants ) {
                log.warn "SNPEFFR: Skipping VCF '${meta.id}' - no variant records found. SNPEFFR requires at least one variant."
                return false
            }
            return true
        }
        .map { id, meta, gz, hasVariants -> [ meta, gz ] }

    ch_tabix_vcf_with_variants.ifEmpty {
        log.warn "SNPEFFR: No VCF files with variant records found. Skipping SNPEFFR analysis."
    }

    SNPEFFR (
        ch_tabix_vcf_with_variants,
        ref_name
    )
    ch_versions = ch_versions.mix(SNPEFFR.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
