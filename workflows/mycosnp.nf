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
if (params.skip_samples_file) { // check for skip_samples_file
    checkPathParamList.add(params.skip_samples_file)
}
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
include { INPUT_CHECK      } from '../subworkflows/local/input_check'
include { BWA_PREPROCESS   } from '../subworkflows/local/bwa-pre-process'
include { BWA_REFERENCE    } from '../subworkflows/local/bwa-reference'
include { GATK_VARIANTS    } from '../subworkflows/local/gatk-variants'
include { CREATE_PHYLOGENY } from '../subworkflows/local/phylogeny'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { GATK4_HAPLOTYPECALLER       } from '../modules/nf-core/modules/gatk4/haplotypecaller/main'
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


workflow MYCOSNP {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

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
    BWA_REFERENCE(ch_fasta)
    ch_versions = ch_versions.mix(BWA_REFERENCE.out.versions)
    fas_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ fa1 ]}.first()
    fai_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ fai ]}.first()
    bai_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ bai ]}.first()
    dict_file = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ dict ]}.first()
    meta_val = BWA_REFERENCE.out.reference_combined.map{meta1, fa1, fai, bai, dict -> [ meta1 ]}.first()

/*
========================================================================================
                          SUBWORKFLOW: Run BWA_PRE_PROCESS
    take:
        reference          channel:  [ tuple reference_fasta, samtools_faidx, bwa_index ]
        reads              channel:  [ val(meta), [ fastq ] ]
    emit:
        alignment           channel: [ val(meta), bam ]
        alignment_index     channel: [ val(meta), bai ]
        alignment_combined  channel: [ val(meta), bam, bai ]
        qualimap            channel: [ val(meta), results ]
        stats               channel: [ val(meta), stats ]
        flagstat            channel: [ val(meta), flagstat ]
        idxstats            channel: [ val(meta), idxstats ]
        versions            channel: [ ch_versions ]    
========================================================================================
*/

    BWA_PREPROCESS( [fas_file, fai_file, bai_file ], INPUT_CHECK.out.reads)
    // do we need to collect sample to perform the qc_report process?
    ch_versions = ch_versions.mix(BWA_PREPROCESS.out.versions)

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
if(! params.skip_vcf)
    {
        GATK4_HAPLOTYPECALLER(  BWA_PREPROCESS.out.alignment_combined.map{meta, bam, bai            -> [ meta, bam, bai, [] ] },
                                fas_file,
                                fai_file,
                                dict_file,
                                [],
                                []
        )
        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
        

        ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{meta, vcf ->[ vcf ]  }.collect()
        ch_vcf_idx = GATK4_HAPLOTYPECALLER.out.tbi.map{meta, idx ->[ idx ]  }.collect()
        GATK4_LOCALCOMBINEGVCFS(
                                    [id:'combined', single_end:false],
                                    ch_vcf,
                                    ch_vcf_idx,
                                    fas_file, 
                                    fai_file, 
                                    dict_file)

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
    }

     CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: Run Pre-FastQC 
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

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
    //ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.post_qc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.stats.map{meta, stats -> [stats]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.flagstat.map{meta, stats -> [stats]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.idxstats.map{meta, stats -> [stats]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.qualimap.map{meta, stats -> [stats]}.ifEmpty([]))
    

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
