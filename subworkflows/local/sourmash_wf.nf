// Sourmash subworkflow

include { SOURMASH_SKETCH             } from '../../modules/local/sourmash_sketch'
include { SOURMASH_COMPARE            } from '../../modules/local/sourmash_compare'
include { MODIFY_MATRIX               } from '../../modules/local/modify_matrix'
include { BEST_MATCH                  } from '../../modules/local/best_match'
include { REPORT_SOURMASH             } from '../../modules/local/report_sourmash'
workflow SOURMASH_WF {

take: ch_spades_contigs

main:
ch_versions=Channel.empty()
//
//Module:
//
SOURMASH_SKETCH(
    ch_spades_contigs
)

ch_signatures = Channel.empty()
ch_signatures = ch_signatures.mix(SOURMASH_SKETCH.out.signatures)
ch_versions = ch_versions.mix(SOURMASH_SKETCH.out.versions)
ch_db = Channel.fromPath('assets/sig_files/*.sig').toList()

//
//Module:
//
SOURMASH_COMPARE(
    ch_signatures,
    ch_db,
    false,
    true
)
ch_sourmash_ani = Channel.empty()
ch_sourmash_ani = ch_sourmash_ani.mix(SOURMASH_COMPARE.out.ani_csv)
ch_versions = ch_versions.mix(SOURMASH_COMPARE.out.versions)
//
//Module:
//
MODIFY_MATRIX(
    ch_sourmash_ani
)
ch_modified_mat = Channel.empty()
ch_modified_mat = ch_modified_mat.mix(MODIFY_MATRIX.out.tsv)
ch_versions = ch_versions.mix(MODIFY_MATRIX.out.versions)
//
//Module:
//
BEST_MATCH(
    ch_modified_mat
)

ch_matches = Channel.empty()
ch_matches = ch_matches.mix(BEST_MATCH.out.match_files.toList())
ch_versions = ch_versions.mix(BEST_MATCH.out.versions)
//
//Module:
//
REPORT_SOURMASH(
    ch_matches
)
ch_report=Channel.empty()
ch_report=ch_report.mix(REPORT_SOURMASH.out.final_report)
ch_versions = ch_versions.mix(REPORT_SOURMASH.out.versions)

emit:
sm_report = ch_report //channel: path(report)
versions  = ch_versions // channel: [ ch_versions ]
}