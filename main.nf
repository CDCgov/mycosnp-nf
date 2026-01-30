#!/usr/bin/env nextflow
/*
========================================================================================
    CDCgov/mycosnp-nf
========================================================================================
    Github : https://github.com/CDCgov/mycosnp-nf
    Wiki   : https://github.com/CDCgov/mycosnp-nf/wiki
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PRE-MYCOSNP
========================================================================================
*/

include { PRE_MYCOSNP_WF          } from './workflows/premycosnp'
include { MYCOSNP                 } from './workflows/mycosnp'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mycosnp_pipeline/main'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mycosnp_pipeline/main'

//
// WORKFLOW: Run pre-mycosnp pipeline
//
workflow PRE_MYCOSNP {
    take:
    samplesheet

    main:
    PRE_MYCOSNP_WF (
        samplesheet
    )

    emit:
    versions = PRE_MYCOSNP_WF.out.versions
    multiqc_report = PRE_MYCOSNP_WF.out.multiqc_report

}

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/


//
// WORKFLOW: Run main CDCgov/mycosnp-nf analysis pipeline
//
workflow NFCORE_MYCOSNP {
    take:
    samplesheet // Define the input parameter

    main:
    MYCOSNP (
        samplesheet
    )

    emit:
    multiqc_report = MYCOSNP.out.multiqc_report
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
    main:
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    // Create a default empty channel for multiqc_report
    multiqc_report = Channel.empty()

    if (params.mode == 'PRE_MYCOSNP') {
        PRE_MYCOSNP(Channel.fromPath(params.input))
        multiqc_report = PRE_MYCOSNP.out.multiqc_report.ifEmpty([])
    } else if (params.mode == 'NFCORE_MYCOSNP') {
        NFCORE_MYCOSNP(Channel.fromPath(params.input))
        multiqc_report = NFCORE_MYCOSNP.out.multiqc_report.ifEmpty([])
    } else {
        log.error "Invalid workflow specified. Use 'PRE_MYCOSNP' or 'NFCORE_MYCOSNP'."
        exit 1
    }

    // Make sure multiqc_report is always valid
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
        multiqc_report.ifEmpty([]) // Ensure it's never null
    )
}

/*
========================================================================================
    THE END
========================================================================================
*/
