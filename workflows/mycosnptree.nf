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
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
include { BWA_REFERENCE    } from '../subworkflows/local/bwa-reference'
include { CREATE_PHYLOGENY } from '../subworkflows/local/phylogeny'
include { GATK_VARIANTS    } from '../subworkflows/local/gatk-variants'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { GATK4_COMBINEGVCFS          } from '../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { SEQKIT_REPLACE              } from '../modules/nf-core/modules/seqkit/replace/main'
include { GATK4_LOCALCOMBINEGVCFS     } from '../modules/local/gatk4_localcombinegvcfs.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []


workflow MYCOSNPTREE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    Channel.fromPath(params.input)
        .splitCsv( header: ['sample', 'vcf', 'vcf_idx'], sep:',', skip:1 )
        .map { row -> [[id:row.sample], row.vcf, row.vcf_idx] }
        .set{ ch_input }

/*
========================================================================================
                          SUBWORKFLOW: Run BWA_REFERENCE
    take:
        fasta               file
    emit:
        masked_fasta        channel: [ val(meta), fas ]
        samtools_index      channel: [ val(meta), fai ]
        bwa_index           channel: [ val(meta), bwa ]
        dict                channel: [ val(meta), dict ]
        reference_combined  channel: [ val(meta), fa, fai, bai, dict ]
        versions            channel: [ ch_versions ]
========================================================================================
*/
/*
    ref_dir                     = null
    ref_fasta                   = null
    ref_fai                     = null
    ref_bwa                     = null
    ref_dict                    = null
*/
    fas_file = Channel.empty()
    fai_file = Channel.empty()
    bai_file = Channel.empty()
    dict_file = Channel.empty()
    meta_val = Channel.empty()


    if(params.ref_dir != null)
    {
        fas_file  = Channel.fromPath(params.ref_dir + "/masked/*.fa*", checkIfExists:true).first()
        fai_file  = Channel.fromPath(params.ref_dir + "/fai/*.fai", checkIfExists:true).first()
        bai_file  = Channel.fromPath(params.ref_dir + "/bwa/bwa", checkIfExists:true, type: 'dir').first()
        dict_file = Channel.fromPath(params.ref_dir + "/dict/*.dict", checkIfExists:true).first()
        // meta_val // Not used
         
    } else if (params.ref_masked_fasta && params.ref_fai && params.ref_bwa && params.ref_dict ) 
    {

        if(params.ref_masked_fasta != null)
        {
            fas_file  = Channel.fromPath(params.ref_masked_fasta, checkIfExists:true).first()
        }
        if(params.ref_fai != null)
        {
            fai_file  = Channel.fromPath(params.ref_fai, checkIfExists:true).first()
        }
        if(params.ref_bwa != null)
        {
            bai_file  = Channel.fromPath(params.ref_bwa, checkIfExists:true, type: 'dir').first()
        }
        if(params.ref_dict != null)
        {
            dict_file = Channel.fromPath(params.ref_dict, checkIfExists:true).first()
        }
    }
    else if (params.fasta) 
    {
        ch_fasta = file(params.fasta) 
        BWA_REFERENCE(ch_fasta)
    
        ch_versions = ch_versions.mix(BWA_REFERENCE.out.versions)
        fas_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ fa1 ]}.first()
        fai_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ fai ]}.first()
        bai_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ bai ]}.first()
        dict_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ dict ]}.first()
        meta_val = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ meta1 ]}.first()
    } else 
    {
        exit 1, 'Input reference fasta or index files not specified!'
    }

/*
========================================================================================
                          SUBWORKFLOW: Run GATK And GATK_VARIANTS
    take:
        fasta
        fai
        bai
        dict
        thismeta
        vcffile
        vcfidx
    emit:
        snps_fasta   channel: [ val(meta), fasta ]
========================================================================================
*/

    ch_vcfs = ch_input.map{meta, vcf, vcf_idx ->[ vcf, vcf_idx ]  }.collect()

    GATK4_LOCALCOMBINEGVCFS(
                                [id:'combined-tree', single_end:false],
                                ch_vcfs,
                                fas_file, 
                                fai_file, 
                                dict_file
                            )

    GATK_VARIANTS( 
                    fas_file, 
                    fai_file, 
                    bai_file, 
                    dict_file,
                    GATK4_LOCALCOMBINEGVCFS.out.combined_gvcf.map{meta, vcf, tbi->[ meta ]}, 
                    GATK4_LOCALCOMBINEGVCFS.out.gvcf, 
                    GATK4_LOCALCOMBINEGVCFS.out.tbi 
                )
    

    ch_versions = ch_versions.mix(GATK_VARIANTS.out.versions)

/*
========================================================================================
                          SUBWORKFLOW: Create Phylogeny 
    take:
        fasta                     file
        constant_sites_string     val: string of constant sites A,C,G,T
    emit:
        rapidnj_tree      = rapidnj_tree     // channel: [ phylogeny ]
        fasttree_tree     = fasttree_tree    // channel: [ phylogeny ]
        iqtree_tree       = iqtree_tree      // channel: [ phylogeny ]
        raxmlng_tree      = raxmlng_tree     // channel: [ phylogeny ]
        versions          = ch_versions 
========================================================================================
*/
    SEQKIT_REPLACE(GATK_VARIANTS.out.snps_fasta) // Swap * for -
    CREATE_PHYLOGENY(SEQKIT_REPLACE.out.fastx.map{meta, fas->[fas]}, '')
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
    
    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
