#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/mycosnp
========================================================================================
    Github : https://github.com/nf-core/mycosnp
    Website: https://nf-co.re/mycosnp
    Slack  : https://nfcore.slack.com/channels/mycosnp
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { MYCOSNP } from './workflows/mycosnp'
include { INPUT_CHECK } from './subworkflows/local/input_check'
include { BWA_PREPROCESS } from './subworkflows/local/bwa-pre-process'
include { BWA_REFERENCE } from './subworkflows/local/bwa-reference'
include { GATK_VARIANTS } from './subworkflows/local/gatk-variants'

//
// WORKFLOW: Run main nf-core/mycosnp analysis pipeline
//
workflow NFCORE_MYCOSNP {
    MYCOSNP ()
}

workflow PREPARE_REFERENCE {
    BWA_REFERENCE ()
}

workflow PREPROCESS {
    BWA_PREPROCESS ()
}

workflow FIND_VARIANTS {
    GATK_VARIANTS ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_MYCOSNP ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
