/*
========================================================================================
    GATK Variants Sub-Workflow
========================================================================================
*/


include { GATK4_HAPLOTYPECALLER } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_COMBINEGVCFS } from '../../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { GATK4_GENOTYPEGVCFS } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTFILTRATION } from '../../modules/nf-core/modules/gatk4/variantfiltration/main'

// TODO: broad vcf filter local module
//include { BROAD_VCFFILTER } from '../../modules/local/broad_vcffilter.nf'

// TODO: broad split vcf local module ? (uses bcftools view, bcftools index, and shell commands)
include { BCFTOOLS_VIEW } from '../../modules/nf-core/modules/bcftools/view/main'
include { BCFTOOLS_QUERY } from '../../modules/nf-core/modules/bcftools/query/main'

include { GATK4_SELECTVARIANTS } from '../../modules/nf-core/modules/gatk4/selectvariants/main'

// TODO : vcf2fasta local module  --> HS: renamed from snp2fasta ==> vcf2fasta
//include { VCFTOFASTA } from '../../modules/local/vcftofasta.nf'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/modules/bcftools/consensus/main'
// TODO: vcf qc report local module
//include { VCF_QCREPORT } from '../../modules/local/vcfqcreport.nf'

combined_gvcf = Channel.empty()

workflow GATK_VARIANTS {

    take:
    reference // channel: [ tuple reference_fasta, fai, bai, dict ]
    combined_gvcf // channel: combined_vcf_file
    

    main:
    ch_versions = Channel.empty()
    

  /*  input = [ [ id:'combined', single_end:false ], 
              combined_gvcf[0], 
              combined_gvcf[1], 
              [], 
              [] 
            ]
  */
    //tuple val(meta), path(gvcf), path(gvcf_index), path(intervals), path(intervals_index)
    //path  fasta
    //path  fasta_index
    //path  fasta_dict
    //path  dbsnp
    //path  dbsnp_index
    GATK4_GENOTYPEGVCFS(
        [ combined_gvcf[0], combined_gvcf[1], combined_gvcf[2], [], [] ],
        reference[0], 
        reference[1], 
        reference[3], [], []
        )
            
    GATK4_VARIANTFILTRATION(
                                 GATK4_GENOTYPEGVCFS.out.vcf.combine(GATK4_GENOTYPEGVCFS.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi]},
                                reference[0], 
                                reference[1], 
                                reference[3]  
                            )

    // gatk  SelectVariants --variant "/data1/vcf-filter.vcf" --reference "/data2/indexed_reference/indexed_reference.fasta" --select-type-to-include "SNP" --output "/data4/gatk-selectvariants/gatk-selectvariants.vcf" 
    GATK4_SELECTVARIANTS(
                                GATK4_VARIANTFILTRATION.out.vcf.combine(GATK4_VARIANTFILTRATION.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi]}
                        )
// Uses
//  docker://geneflow/python:2.7.18-scipy python
// Uses library vcftools in assets directory
  // filterGatkGenotypes.py --min_GQ "50" 
  //                        --keep_GQ_0_refs 
  //                        --min_percent_alt_in_AD "0.8" 
  //                        --min_total_DP "10" 
  //                        --keep_all_ref "/data7/gatk-variantfiltration.vcf" 
  //                         > "vcf-filter.vcf"
    // TODO // BROAD_VCF_FILTER()



    // TODO //SPLITVCF() split-vcf-selectvariants --> also uses bcftools view, bcftools index, and shell commands + the GATK4 output
    //BCFTOOLS_VIEW()
    //BCFTOOLS_QUERY()


    // TODO //VCFTOFASTA()


    //BCFTOOLS_CONSENSUS()

    // TODO //VCF_QCREPORT()


// Collect versions information

    ch_versions = ch_versions.mix(  GATK4_GENOTYPEGVCFS.out.versions, 
                                    GATK4_VARIANTFILTRATION.out.versions, 
                                    GATK4_SELECTVARIANTS.out.versions
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
