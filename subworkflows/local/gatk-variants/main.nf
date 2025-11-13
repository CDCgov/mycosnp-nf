/*
========================================================================================
    GATK Variants Sub-Workflow
========================================================================================
*/


include { GATK4_GENOTYPEGVCFS                     } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTFILTRATION                 } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_SELECTVARIANTS                    } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { FILTER_GATK_GENOTYPES                   } from '../../../modules/local/vcftools/main'
include { BCFTOOLS_INDEX                          } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_VIEW  as BCFTOOLS_VIEW_CONVERT } from '../../../modules/nf-core/bcftools/view/main'
include { SPLIT_VCF                               } from '../../../modules/local/splitvcf/main'
include { VCF_TO_FASTA                            } from '../../../modules/local/vcftofasta/main'
include { VCF_QC                                  } from '../../../modules/local/vcfqc/main'
include { BCFTOOLS_QUERY                          } from '../../../modules/nf-core/bcftools/query/main'
include { VCF_CONSENSUS                           } from '../../../modules/local/vcfconsensus/main'

workflow GATK_VARIANTS {

    take:
    fasta               // channel: tuple val(meta), path(fasta)
    fai                 // channel: tuple val(meta), path(fai)
    bai                 // channel: tuple val(meta), path(bai)
    dict                // channel: tuple val(meta), path(dict)
    combined_vcfidx     // channels: tuple val(meta), path(vcf), path(vcfidx)
    max_amb_samples     // val: maximum number of ambiguous samples allowed per variant
                            // to be retained in final VCF
                            // (e.g. if 5 samples, and max_amb_samples=2, variants with >2 amb samples are filtered out)
                            // set to -1 to disable filtering based on ambiguous samples
    max_perc_amb_samples // val: maximum percentage of ambiguous samples allowed per variant
                             // to be retained in final VCF
                             // (e.g. if 10 samples, and max_perc_amb_samples=20, variants with >2 amb samples are filtered out)
                             // set to 0 to disable filtering based on percentage of ambiguous samples
    min_depth            // val: minimum depth required to call a SNP at a given position

    main:
    ch_versions = Channel.empty()

    GATK4_GENOTYPEGVCFS (
        combined_vcfidx.map { meta, vcf, idx -> [ meta, vcf, idx, [], [] ] },
        fasta.map { f -> [[:], f] },
        fai.map { f -> [[:], f] },
        dict.map { d -> [[:], d] },
        Channel.value([[:], []]),
        Channel.value([[:], []])
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)


    GATK4_VARIANTFILTRATION (
        GATK4_GENOTYPEGVCFS.out.vcf.combine(GATK4_GENOTYPEGVCFS.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi]},
        fasta.map { f -> [[:], f] },
        fai.map { f -> [[:], f] },
        dict.map { d -> [[:], d] },
        Channel.value([[:], []])
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions)

    GATK4_SELECTVARIANTS (
        GATK4_VARIANTFILTRATION.out.vcf.combine(GATK4_VARIANTFILTRATION.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi, []]}
    )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)

    FILTER_GATK_GENOTYPES (
        GATK4_SELECTVARIANTS.out.vcf
    )
    ch_versions = ch_versions.mix(FILTER_GATK_GENOTYPES.out.versions)

    fin_comb_vcf = FILTER_GATK_GENOTYPES.out.vcf.first()

    BCFTOOLS_VIEW_CONVERT (
        fin_comb_vcf.map { meta, vcf->[ meta, vcf, [] ] }, [], [], []
    )   // Convert to bgzip
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_CONVERT.out.versions)

    BCFTOOLS_INDEX (
        BCFTOOLS_VIEW_CONVERT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    SPLIT_VCF (
        BCFTOOLS_VIEW_CONVERT.out.vcf.combine(BCFTOOLS_INDEX.out.csi).map{meta1, vcf, meta2, csi-> [meta1, vcf, csi] }
    )
    ch_versions = ch_versions.mix(SPLIT_VCF.out.versions)

    fin_comb_vcf.combine (
        SPLIT_VCF.out.txt
        ).map {
            meta1, vcf, meta2, txt ->
            [ meta1, vcf, txt, max_amb_samples, max_perc_amb_samples, min_depth ]
        }.set { final_vcf_txt }

    VCF_CONSENSUS (
        BCFTOOLS_VIEW_CONVERT.out.vcf.combine ( BCFTOOLS_INDEX.out.csi ).map { meta1, vcf, meta2, csi-> [meta1, vcf, csi] },
        fasta
    )
    ch_versions = ch_versions.mix(VCF_CONSENSUS.out.versions)

    VCF_TO_FASTA (
        final_vcf_txt,
        fasta
    )
    ch_versions = ch_versions.mix(VCF_TO_FASTA.out.versions)

    VCF_QC (
        VCF_TO_FASTA.out.fasta.map { meta, fa-> fa }
    )
    ch_versions = ch_versions.mix(VCF_QC.out.versions)

    emit:
    snps_fasta   = VCF_TO_FASTA.out.fasta        // channel: [ val(meta), fasta ]
    versions     = ch_versions                   // channel: [ versions.yml ]
    filtered_vcf = BCFTOOLS_VIEW_CONVERT.out.vcf // channel: [ val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}") ]
    // filtered_vcf    = BROAD_VCFFILTER.out
    // split_vcf_broad = SPLITVCF.out        --> the broad vcf file
    // variants        = GATK4_SELECTVARIANTS.out
    // split_vcf_gatk4 = SPLITVCF.out        --> the gatk4 vcf file
    // consensus_fasta = BCFTOOLS_CONSENSUS.out
    // qc_report       = VCF_QCREPORT.out

}
