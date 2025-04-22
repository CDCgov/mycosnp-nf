#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/mycosnp
========================================================================================
    Github : https://github.com/CDCgov/mycosnp-nf
    Wiki   : https://github.com/CDCgov/mycosnp-nf/wiki
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
    PRE-MYCOSNP
========================================================================================
*/

include { PRE_MYCOSNP_WF } from './workflows/pre_mycosnp'

//
// WORKFLOW: Run pre-mycosnp pipeline
//
workflow PRE_MYCOSNP {
    PRE_MYCOSNP_WF ()
}

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { MYCOSNP } from './workflows/mycosnp'

//
// WORKFLOW: Run main nf-core/mycosnp analysis pipeline
//
workflow NFCORE_MYCOSNP {
    MYCOSNP ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute the specified workflow
//
workflow {
    if (params.workflow == 'PRE_MYCOSNP') {
        PRE_MYCOSNP()
    } else if (params.workflow == 'NFCORE_MYCOSNP') {
        NFCORE_MYCOSNP()
    } else {
        log.error "Invalid workflow specified. Use 'PRE_MYCOSNP' or 'NFCORE_MYCOSNP'."
        exit 1
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
