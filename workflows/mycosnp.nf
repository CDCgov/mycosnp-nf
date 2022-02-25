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

    // SUBWORKFLOW: Run BWA_PRE_PROCESS
    // take:
    // tuple reference_fasta, samtools_faidx, bwa_index
    // tuple meta, fastq
    
    BWA_PREPROCESS( [BWA_REFERENCE.out.masked_fasta, BWA_REFERENCE.out.samtools_index, BWA_REFERENCE.out.bwa_index ], INPUT_CHECK.out.reads)
    //BWA_PRE_PROCESS(BWA_REFERENCE.out.map{masked_fasta, samtools_index, bwa_index, dict, versions->[masked_fasta, samtools_index, bwa_index] }, INPUT_CHECK.out.reads)
    //ch_versions = ch_versions.mix(BWA_PRE_PROCESS.out.versions)

    // SUBWORKFLOW: Run GATK_VARIANTS
    // tuple reference_fasta, samtools_faidx, bwa_index
    // tuple meta, alignment, aligment_index
    // emit:
    // alignment = ADDORREPLACEGROUPS.out  -> (meta, bam)
    // alignment_index = BAM_INDEX.out (meta, bam_index)
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

    //ch_gatk_in = ch_gatk_in.mix(BWA_PREPROCESS.out.alignment)
    //ch_gatk_in = ch_gatk_in.mix(BWA_PREPROCESS.out.alignment_index)
    //ch_grouped = Channel.empty()
    //ch_grouped = ch_gatk_in.groupTuple()  
    //ch_gatk_in.view()
    //ch_grouped.view()
    //ch_aln = Channel.empty()
    //ch_aln = BWA_PREPROCESS.out.alignment.collect()
    //ch_aln_idx = Channel.empty()
    //ch_aln_idx = BWA_PREPROCESS.out.alignment_index.collect()
    //ch_mixed = Channel.empty()
    
    //ch_mixed = ch_mixed.mix(ch_aln)
    //ch_mixed = ch_mixed.mix(ch_aln_idx)
    
    //ch_grouped = Channel.empty()
    //ch_grouped = ch_mixed.groupTuple()
    //ch_mixed.view()
    // ch_grouped.view()

    ch_all_aln = Channel.empty()
    ch_all_aln = BWA_PREPROCESS.out.alignment_combined.collect()
    ch_all_aln.view()

    GATK_VARIANTS( [BWA_REFERENCE.out.masked_fasta, BWA_REFERENCE.out.samtools_index, BWA_REFERENCE.out.bwa_index ], ch_all_aln )
    //GATK_VARIANTS( [BWA_REFERENCE.out.masked_fasta, BWA_REFERENCE.out.samtools_index, BWA_REFERENCE.out.bwa_index ], ch_gatk_in )

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
