/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config
    // params.snpeffdb
]
// check for skip_samples_file
if (params.skip_samples_file) { checkPathParamList.add(params.skip_samples_file) }

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.input) { samplesheet = Channel.fromPath(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet file not specified!' }


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
include { SRA_FASTQ_SRATOOLS       } from '../subworkflows/local/sra_fastq_sratools/main'
include { INPUT_CHECK              } from '../subworkflows/local/input_check/main'
include { SEQKIT_PAIR              } from '../modules/nf-core/seqkit/pair/main'
include { FAQCS                    } from '../modules/local/faqcs/main'
include { GAMBIT_QUERY             } from '../modules/local/gambit/main'
include { SUBTYPE                  } from '../modules/local/subtype/main'
include { PRE_MYCOSNP_INDV_SUMMARY } from '../modules/local/pre_mycosnp_indv_summary/main'
include { PRE_MYCOSNP_COMB_SUMMARY } from '../modules/local/pre_mycosnp_comb_summary/main'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/fastqc/main'
include { SHOVILL as SHOVILL          } from '../modules/local/shovill/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary


workflow PRE_MYCOSNP_WF {

    take:
    samplesheet                 // New samplesheet combines ingestion for fastq reads, sra accessions, and vcf files for phylogeny

    main:

    def ch_versions = Channel.empty()

    // Create empty channels for reference files required by the main workflow
    def fas_file  = Channel.empty()
    def fai_file  = Channel.empty()
    def bai_file  = Channel.empty()
    def dict_file = Channel.empty()
    def meta_val  = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    def ch_all_reads = Channel.empty()
    def ch_sra_list  = Channel.empty()
    INPUT_CHECK( params.input )
    ch_all_reads = ch_all_reads.mix(INPUT_CHECK.out.ch_fastq_reads)       // channel: [ val(meta), [ reads ] ]
    ch_sra_list  = ch_sra_list.mix (INPUT_CHECK.out.ch_sra_list)          // channel: [ val(meta), sra_id    ]
    ch_versions  = ch_versions.mix (INPUT_CHECK.out.versions)             // channel: [ versions.yml         ]
    // Other output channels from INPUT_CHECK for VCFs are not used in this wf

    //
    // SUBWORKFLOW: Fetch FASTQ reads from input SRA accession IDs, mix with fastq reads from samplesheet
    //
    SRA_FASTQ_SRATOOLS(ch_sra_list)
    ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)
    ch_versions  = ch_versions.mix(SRA_FASTQ_SRATOOLS.out.versions)

    //
    // MODULE: Run Pre-FastQC
    //
    FASTQC_RAW (
        ch_all_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    //
    // MODULE: Run seqkit to remove unpaired reads
    //
    SEQKIT_PAIR(
        ch_all_reads
    )
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions.first())

    //
    // MODULE: Run FAQCs - no downsampling option because a reference cannot be supplied before knowing the species
    //
    FAQCS(
        SEQKIT_PAIR.out.reads
    )
    ch_versions = ch_versions.mix(FAQCS.out.versions.first())

    //
    // MODULE: Run Shovill
    //
    SHOVILL (
        FAQCS.out.reads
    )
    ch_versions = ch_versions.mix(SHOVILL.out.versions.first())

    //
    // MODULE: Run Gambit
    //
    GAMBIT_QUERY(
        SHOVILL.out.contigs,
        params.gambit_db,
        params.gambit_h5_dir
    )
    ch_versions = ch_versions.mix(GAMBIT_QUERY.out.versions.first())

    // Join the GAMBIT output and the spades assembly into a single channel
    SHOVILL.out.contigs.map  { meta, contigs -> [meta, contigs] }.set{ ch_contigs }
    GAMBIT_QUERY.out.taxa.map{ meta, gambit  -> [meta, gambit]  }
        .join(ch_contigs)
        .set{ ch_gambit_assembly }

    //
    // MODULE: Subtype
    //
    SUBTYPE(
        ch_gambit_assembly,
        params.subtype_db
    )
    ch_versions = ch_versions.mix(SUBTYPE.out.versions.first())

    //
    // TODO? --> MODULE: Create line summary for each sample
    //

    // Combine trimmed reads and the QC reference into single channel
    FAQCS.out.txt.map        { meta, txt     -> [meta, txt    ] }.set{ ch_faqcs_txt }
    GAMBIT_QUERY.out.taxa.map{ meta, gambit  -> [meta, gambit ] }.set{ ch_gambit    }
    SUBTYPE.out.subtype.map  { meta, subtype -> [meta, subtype] }.set{ ch_subtype   }
    SHOVILL.out.contigs.map  { meta, contigs -> [meta, contigs] }
        .join(ch_faqcs_txt)
        .join(ch_gambit)
        .join(ch_subtype)
        .set{ ch_line_summary_input }

    PRE_MYCOSNP_INDV_SUMMARY(
        ch_line_summary_input
    )
    ch_versions = ch_versions.mix(PRE_MYCOSNP_INDV_SUMMARY.out.versions)
    //
    // MODULE: Combine line summaries into single output
    //
    PRE_MYCOSNP_COMB_SUMMARY(
        PRE_MYCOSNP_INDV_SUMMARY.out.result
            .map{ meta, result -> [result] }
            .collect()
    )
    ch_versions = ch_versions.mix(PRE_MYCOSNP_COMB_SUMMARY.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique()
        .collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        [],
        [],
        []
    )
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    emit:
    multiqc_report = MULTIQC.out.report
    versions = ch_versions
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/
/*
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}
*/
/*
========================================================================================
    THE END
========================================================================================
*/
