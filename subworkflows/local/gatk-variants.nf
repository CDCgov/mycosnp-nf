/*
========================================================================================
    GATK Variants Sub-Workflow
========================================================================================
*/


include { GATK4_HAPLOTYPECALLER } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_COMBINEGVCFS } from '../../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { GATK4_GENOTYPEGVCFS } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTFILTRATION } from '../../modules/nf-core/modules/gatk4/variantfiltration/main'
include { GATK4_SELECTVARIANTS } from '../../modules/nf-core/modules/gatk4/selectvariants/main'

// TODO: broad vcf filter local module
include { FILTER_GATK_GENOTYPES } from '../../modules/local/vcftools.nf'

// TODO: broad split vcf local module ? (uses bcftools view, bcftools index, and shell commands)
include { BCFTOOLS_INDEX } from '../../modules/nf-core/modules/bcftools/index/main'
include { BCFTOOLS_VIEW  as BCFTOOLS_VIEW_CONVERT } from '../../modules/nf-core/modules/bcftools/view/main'
include { SPLIT_VCF } from '../../modules/local/splitvcf.nf'
include { VCF_TO_FASTA } from '../../modules/local/vcftofasta.nf'

include { BCFTOOLS_QUERY } from '../../modules/nf-core/modules/bcftools/query/main'
include { VCF_CONSENSUS } from '../../modules/local/vcfconsensus.nf'
//include { VCF_QCREPORT } from '../../modules/local/vcfqcreport.nf'

combined_gvcf = Channel.empty()

workflow GATK_VARIANTS {

    take:
    //reference // channel: [ tuple reference_fasta, fai, bai, dict ]
    fasta
    fai
    bai
    dict
    thismeta
    vcffile
    vcfidx

    main:
    ch_versions = Channel.empty()
    //combined_gvcf = Channel.empty()
    combined_gvcf = thismeta.combine(vcffile).combine(vcfidx)
    combined_gvcf.view()

   
    GATK4_GENOTYPEGVCFS(
        combined_gvcf.map{meta, vcf, idx -> [ meta, vcf, idx, [], [] ]},
        fasta, 
        fai, 
        dict, [], []
        )
       // GATK4_GENOTYPEGVCFS.out.vcf.view()

    GATK4_VARIANTFILTRATION(
                                 GATK4_GENOTYPEGVCFS.out.vcf.combine(GATK4_GENOTYPEGVCFS.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi]},
                                fasta, 
                                fai, 
                                dict
                            )

    GATK4_SELECTVARIANTS(
                                GATK4_VARIANTFILTRATION.out.vcf.combine(GATK4_VARIANTFILTRATION.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi]}
                        )

    FILTER_GATK_GENOTYPES(GATK4_SELECTVARIANTS.out.vcf)
    fin_comb_vcf = FILTER_GATK_GENOTYPES.out.vcf.first()

    // Convert to bgzip
    BCFTOOLS_VIEW_CONVERT(fin_comb_vcf.map{meta, vcf->[ meta, vcf, [] ] }, [], [], []  )
    BCFTOOLS_INDEX(BCFTOOLS_VIEW_CONVERT.out.vcf)
    SPLIT_VCF(
                     BCFTOOLS_VIEW_CONVERT.out.vcf.combine(BCFTOOLS_INDEX.out.csi).map{meta1, vcf, meta2, csi-> [meta1, vcf, csi] }
             )


    final_vcf_txt = Channel.empty()
    fin_comb_vcf.combine(SPLIT_VCF.out.txt).map{meta1, vcf, meta2, txt -> 
                    [ meta1, vcf, txt, params.max_amb_samples, params.max_perc_amb_samples]}.set{final_vcf_txt}

    VCF_CONSENSUS(
        BCFTOOLS_VIEW_CONVERT.out.vcf.combine(BCFTOOLS_INDEX.out.csi).map{meta1, vcf, meta2, csi-> [meta1, vcf, csi] },
        fasta
    )

    VCF_TO_FASTA(final_vcf_txt, fasta)


    // TODO //VCF_QCREPORT()


// Collect versions information

    ch_versions = ch_versions.mix(  GATK4_GENOTYPEGVCFS.out.versions, 
                                    GATK4_VARIANTFILTRATION.out.versions, 
                                    GATK4_SELECTVARIANTS.out.versions,
                                    BCFTOOLS_VIEW_CONVERT.out.versions,
                                    BCFTOOLS_INDEX.out.versions,
                                    SPLIT_VCF.out.versions,
                                    VCF_CONSENSUS.out.versions
                                )
                                
 

    emit:
    // filtered_vcf = BROAD_VCFFILTER.out
    // split_vcf_broad = SPLITVCF.out        --> the broad vcf file
    // variants = GATK4_SELECTVARIANTS.out
    // split_vcf_gatk4 = SPLITVCF.out        --> the gatk4 vcf file
    // snps_fasta = VCFTOFASTA.out
    // consensus_fasta = BCFTOOLS_CONSENSUS.out
    // qc_report = VCF_QCREPORT.out
    versions = ch_versions // channel: [ versions.yml ]
}

/*

Workflow process:

1. Call variants using the GATK 4.1.4.1 HaplotypeCaller tool.
2. Combine gVCF files from the HaplotypeCaller into a single VCF using the GATK 4.1.4.1 CombineGVCFs tool.
3. Call genotypes using the GATK 4.1.4.1 GenotypeGVCFs tool.
4. Filter the variants using the GATK 4.1.4.1 VariantFiltration tool and the default (but customizable) filter: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10'.
5. Run a customized VCF filtering script provided by the Broad Institute.
6. Split the filtered VCF file by sample.
7. Select only SNPs from the VCF files using the GATK 4.1.4.1 SelectVariants tool.
8. Split the VCF file with SNPs by sample.
9. Create a consensus sequence for each sample using BCFTools 1.9 and SeqTK 1.2.
10. Create a multi-fasta file from the VCF SNP positions using a custom script from Broad.

*/
