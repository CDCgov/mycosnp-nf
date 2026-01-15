/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check input path parameters to see if they exist
/*
def checkPathParamList = [
    params.input,
    params.multiqc_config
]
*/
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SRA_FASTQ_SRATOOLS          } from '../subworkflows/local/sra_fastq_sratools/main'
include { INPUT_CHECK                 } from '../subworkflows/local/input_check/main'
include { GAMBIT_QUERY                } from '../modules/local/gambit/main'
include { SUBTYPE                     } from '../modules/local/subtype/main'
include { PRE_MYCOSNP_INDV_SUMMARY    } from '../modules/local/pre_mycosnp_indv_summary/main'
include { SHOVILL as SHOVILL          } from '../modules/local/shovill/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SEQKIT_PAIR                 } from '../modules/nf-core/seqkit/pair/main'
include { FAQCS                       } from '../modules/nf-core/faqcs/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow PRE_MYCOSNP_WF {

    take:
    samplesheet // New samplesheet combines ingestion for fastq reads, sra accessions, and vcf files for phylogeny

    main:
    def ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        samplesheet
    )
    ch_versions  = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_all_reads = INPUT_CHECK.out.ch_fastq_reads // channel: [ val(meta), [ reads ] ]
    //
    // SUBWORKFLOW: Fetch FASTQ reads from input SRA accession IDs, mix with fastq reads from samplesheet
    //
    SRA_FASTQ_SRATOOLS (
        INPUT_CHECK.out.ch_sra_list
    )
    ch_versions  = ch_versions.mix(SRA_FASTQ_SRATOOLS.out.versions)

    ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)
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
    SEQKIT_PAIR (
        ch_all_reads
    )
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions.first())

    //
    // MODULE: Run FAQCs - no downsampling option because a reference cannot be supplied before knowing the species
    //
    FAQCS (
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
    GAMBIT_QUERY (
        SHOVILL.out.contigs,
        params.gambit_db,
        params.gambit_h5_dir
    )
    ch_versions = ch_versions.mix(GAMBIT_QUERY.out.versions.first())

    // Join the GAMBIT output and the spades assembly into a single channel
    SHOVILL.out.contigs.map { meta, contigs -> [meta, contigs] }.set { ch_contigs }
    GAMBIT_QUERY.out.taxa.map { meta, gambit  -> [meta, gambit] }
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

    // Combine trimmed reads and the QC reference into single channel
    // PRE_MYCOSNP_INDV_SUMMARY needs both the stats.txt file and the debug directory from FAQCS
    ch_faqcs_combined = FAQCS.out.stats
        .join(FAQCS.out.debug)

    GAMBIT_QUERY.out.taxa.map{ meta, gambit  -> [meta, gambit ] }.set{ ch_gambit    }
    SUBTYPE.out.subtype.map  { meta, subtype -> [meta, subtype] }.set{ ch_subtype   }
    SHOVILL.out.contigs.map  { meta, contigs -> [meta, contigs] }
        .join(ch_faqcs_combined)
        .join(ch_gambit)
        .join(ch_subtype)
        .set{ ch_line_summary_input }

    PRE_MYCOSNP_INDV_SUMMARY (
        ch_line_summary_input
    )
    ch_versions = ch_versions.mix(PRE_MYCOSNP_INDV_SUMMARY.out.versions)

    //
    // Combine line summaries into single output using collectFile
    //
    ch_combined_summary = PRE_MYCOSNP_INDV_SUMMARY.out.result
        .map { meta, result -> result }
        .collectFile(
            name: 'pre-mycosnp-summary.csv',
            storeDir: "${params.outdir}/aggregate_outputs",
            seed: "sample,PM_Predicted_Rank,PM_Predicted_Taxon,PM_Subtype_Closest_Match,PM_Subtype_ANI,PM_Closest_GAMBIT_Entry_Description,PM_Closest_GAMBIT_Entry_Distance,PM_Trimmed_Reads,PM_Avg_Read_Quality,PM_Sample_Assembly_Length,PM_Sample_Assembly_GC,PM_Reference_Genome_Length,PM_Avg_Depth_Coverage,PM_Reference_GC\n",
            newLine: false,
            sort: true
        )

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'premycosnp_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //

    ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.map{it[1]})

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_custom_config.toList(),
        [],
        [],
        []
    )
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    emit:
    multiqc_report = MULTIQC.out.report
    versions       = ch_versions
}

/*
========================================================================================
    THE END
========================================================================================
*/
