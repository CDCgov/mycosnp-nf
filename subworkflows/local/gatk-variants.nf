/*
========================================================================================
    GATK Variants Sub-Workflow
========================================================================================
*/


include { GATK4_HAPLOTYPECALLER } from '../../modules/nf-core/modules/gatk4/haplotypecaller/'
include { GATK4_COMBINEGVCFS } from '../../modules/nf-core/modules/gatk4/combinegvcfs/'
include { GATK4_GENOTYPEGVCFS } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/'
include { GATK4_VARIANTFILTRATION } from '../../modules/nf-core/modules/gatk4/variantfiltration/'

// TODO: broad vcf filter local module
//include { BROAD_VCFFILTER } from '../../modules/local/broad_vcffilter.nf'

// TODO: broad split vcf local module ? (uses bcftools view, bcftools index, and shell commands)
include { BCFTOOLS_VIEW } from '../../modules/nf-core/modules/bcftools/view/'
include { BCFTOOLS_QUERY } from '../../modules/nf-core/modules/bcftools/query/'

include { GATK4_SELECTVARIANTS } from '../../modules/nf-core/modules/gatk4/selectvariants/'

// TODO : vcf2fasta local module  --> HS: renamed from snp2fasta ==> vcf2fasta
//include { VCFTOFASTA } from '..../../modules/local/vcftofasta.nf'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/modules/bcftools/consensus/'
// TODO: vcf qc report local module
//include { VCF_QCREPORT } from '..../../modules/local/vcfqcreport.nf'


workflow GATK_VARIANTS {

    take:
    tuple reference_fasta, samtools_faidx, bwa_index
    tuple meta, alignment, aligment_index

    main:
    GATK4_HAPLOTYPECALLER()
    GATK4_COMBINEGVCFS()
    GATK4_GENOTYPEGVCFS()
    GATK4_VARIANTFILTRATION()
    //BROAD_VCFFILTER()
    //SPLITVCF() split-vcf-broad --> uses bcftools view, bcftools index, and shell commands
    BCFTOOLS_VIEW()
    BCFTOOLS_QUERY()
    //GATK4_SELECTVARIANTS()
    //SPLITVCF() split-vcf-selectvariants --> also uses bcftools view, bcftools index, and shell commands + the GATK4 output
    BCFTOOLS_VIEW()
    BCFTOOLS_QUERY()
    //VCFTOFASTA()
    BCFTOOLS_CONSENSUS()
    //VCF_QCREPORT()

    emit:
    //filtered_vcf = BROAD_VCFFILTER.out
    //split_vcf_broad = SPLITVCF.out        --> the broad vcf file
    variants = GATK4_SELECTVARIANTS.out
    //split_vcf_gatk4 = SPLITVCF.out        --> the gatk4 vcf file
    //snps_fasta = VCFTOFASTA.out
    consensus_fasta = BCFTOOLS_CONSENSUS.out
    //qc_report = VCF_QCREPORT.out
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

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
