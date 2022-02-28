/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMycosnp.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Input reference fasta not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { BWA_PREPROCESS } from '../subworkflows/local/bwa-pre-process'
include { BWA_REFERENCE } from '../subworkflows/local/bwa-reference'
include { GATK_VARIANTS } from '../subworkflows/local/gatk-variants'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
//include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { GATK4_HAPLOTYPECALLER }       from '../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_COMBINEGVCFS }          from '../modules/nf-core/modules/gatk4/combinegvcfs/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MYCOSNP {

    ch_versions = Channel.empty()
    ch_gatk_in = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // SUBWORKFLOW: Run BWA_REFERENCE
    //emit:
    //fasta = masked_fasta
    //samtools_index = SAMTOOLS_FAIDX.out.fai
    //bwa_index = BWA_INDEX.out.index
    //dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
    //versions = ch_versions // channel: [ versions.yml ]

    BWA_REFERENCE(ch_fasta)
    ch_versions = ch_versions.mix(BWA_REFERENCE.out.versions)
    fas_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ fa1 ]}.first()
    fai_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ fai ]}.first()
    bai_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ bai ]}.first()
    dict_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ dict ]}.first()
    meta_val = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ meta1 ]}.first()

    // SUBWORKFLOW: Run BWA_PRE_PROCESS
    // take:
    // tuple reference_fasta, samtools_faidx, bwa_index
    // tuple meta, fastq
    
    BWA_PREPROCESS( [fas_file, fai_file, bai_file ], INPUT_CHECK.out.reads)
    ch_versions = ch_versions.mix(BWA_PREPROCESS.out.versions)
    //ch_versions = ch_versions.mix(BWA_PRE_PROCESS.out.versions)

    // SUBWORKFLOW: Run GATK And GATK_VARIANTS

    GATK4_HAPLOTYPECALLER(  BWA_PREPROCESS.out.alignment_combined.map{meta, bam, bai            -> [ meta, bam, bai, [] ] },
                            fas_file,
                            fai_file,
                            dict_file,
                            [],
                            []
     )
     ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
     
  //  ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{meta, vcf ->[ vcf ]  }.collect()
  //  ch_vcf_idx = GATK4_HAPLOTYPECALLER.out.tbi.map{meta, idx ->[ idx ]  }.collect()

  /*
    GATK4_COMBINEGVCFS( tuple( [ id:'test', single_end:false ], 
                          GATK4_HAPLOTYPECALLER.out.vcf.map{meta, vcf ->[ vcf ]  }, 
                          GATK4_HAPLOTYPECALLER.out.tbi.map{meta, idx ->[ idx ]  } )
                        , 
                        fas_file, 
                        fai_file, 
                        dict_file )
    */
    //GATK4_COMBINEGVCFS([ [ id:'test', single_end:false ], ch_vcf, ch_vcf_idx] , fas_file, fai_file, dict_file )
    //GATK4_COMBINEGVCFS(inputvcf , fas_file, fai_file, dict_file )
    //GATK_VARIANTS( [fas_file, fai_file, bai_file, dict_file ], GATK4_COMBINEGVCFS.out.combined_gvcf )



    // These files are temporary until combinegvcfs working - this will allow pipeline to continue with testdata running only
    gvcftest = file("$projectDir/assets/testdata/combinegvcfs/gatk-combinegvcfs.g.vcf")
    gvcfidxtest = file("$projectDir/assets/testdata/combinegvcfs/gatk-combinegvcfs.g.vcf.idx")
    GATK_VARIANTS( [fas_file, fai_file, bai_file, dict_file ], [ [ id:'combined', single_end:false ], gvcftest, gvcfidxtest] )
    ch_versions = ch_versions.mix(GATK_VARIANTS.out.versions)



     CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /*
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMycosnp.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    */
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
