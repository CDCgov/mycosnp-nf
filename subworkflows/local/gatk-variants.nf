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


workflow GATK_VARIANTS {

    take:
    reference // channel: [ tuple meta, reference_fasta, fai, bai, dict ]
    alignments // channel: [ tuple meta, alignment, aligment_index ]
    

    main:
    ch_versions = Channel.empty()

    ch_vcf = Channel.empty()
    //alignments.collate(3, false)
    //ch_vcf = alignments.collate(3)
    //ch_vcf.view()

    // https://gitlab.com/geneflow/apps/gatk-haplotypecaller-gf2.git
    // gatk  HaplotypeCaller --input "/data1/SRR13710812/SRR13710812.bam" 
    //    --sample-ploidy "1" 
    //    --emit-ref-confidence "GVCF" 
    //    --native-pair-hmm-threads "4" 
    //    --reference "/data5/indexed_reference/indexed_reference.fasta" 
    //    --output "/data7/SRR13710812/SRR13710812.g.vcf
    /*
        tuple val(meta), path(input), path(input_index), path(intervals)
    path fasta
    path fai
    path dict
    path dbsnp
    path dbsnp_tbi
    */
    //GATK4_HAPLOTYPECALLER()
    //GATK4_COMBINEGVCFS()
    
    //GATK4_GENOTYPEGVCFS()
    //GATK4_VARIANTFILTRATION()
    // TODO //BROAD_VCFFILTER() 
    // TODO //SPLITVCF() split-vcf-broad --> uses bcftools view, bcftools index, and shell commands
    //BCFTOOLS_VIEW()
    //BCFTOOLS_QUERY()

    // gatk  SelectVariants --variant "/data1/vcf-filter.vcf" --reference "/data2/indexed_reference/indexed_reference.fasta" --select-type-to-include "SNP" --output "/data4/gatk-selectvariants/gatk-selectvariants.vcf" 
    //GATK4_SELECTVARIANTS()

    // TODO //SPLITVCF() split-vcf-selectvariants --> also uses bcftools view, bcftools index, and shell commands + the GATK4 output
    //BCFTOOLS_VIEW()
    //BCFTOOLS_QUERY()
    // TODO //VCFTOFASTA()
    //BCFTOOLS_CONSENSUS()
    // TODO //VCF_QCREPORT()

// Collect versions information
/*
    ch_versions = ch_versions.mix(  NUCMER.out.versions, 
                                    BWA_INDEX.out.versions, 
                                    SAMTOOLS_FAIDX.out.versions, 
                                    PICARD_CREATESEQUENCEDICTIONARY.out.versions
                                )
*/                                
                                
    emit:
    //filtered_vcf = BROAD_VCFFILTER.out
    //split_vcf_broad = SPLITVCF.out        --> the broad vcf file
    // variants = GATK4_SELECTVARIANTS.out
    //split_vcf_gatk4 = SPLITVCF.out        --> the gatk4 vcf file
    //snps_fasta = VCFTOFASTA.out
    // consensus_fasta = BCFTOOLS_CONSENSUS.out
    //qc_report = VCF_QCREPORT.out
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
