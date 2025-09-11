//
// Run Snpeff Build
//

include { SNPEFF_BUILD        } from '../../../modules/local/snpeff_build/main'


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

    emit:
    snpeff_db        = ch_snpeff_db                    // channel: [ val(meta), [ db ] ]
    config           = ch_snpeff_config                // channel: [ val(meta), [ config ] ]
    versions         = ch_versions                     // channel: [ versions.yml ]

}
