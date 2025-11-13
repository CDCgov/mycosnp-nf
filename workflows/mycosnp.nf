/*
========================================================================================
    VALIDATE MYCOSNP WORKFLOW SPECIFIC INPUTS
========================================================================================
*/

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta
    // params.snpeffdb
]

// check for skip_samples_file
if (params.skip_samples_file) { checkPathParamList.add(params.skip_samples_file) }

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet file not specified!' }

// if (params.snpeffdb == null) { exit 1, 'Input path to snpeffdb not specified!' }
// if (params.snpeffconfig == null) { exit 1, 'Input snpeff config file not specified' }


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

include { SRA_FASTQ_SRATOOLS   } from '../subworkflows/local/sra_fastq_sratools/main'
include { INPUT_CHECK          } from '../subworkflows/local/input_check/main'
include { BWA_PREPROCESS       } from '../subworkflows/local/bwa-pre-process/main'
include { BWA_REFERENCE        } from '../subworkflows/local/bwa-reference/main'
include { GATK_VARIANTS        } from '../subworkflows/local/gatk-variants/main'
include { CREATE_PHYLOGENY     } from '../subworkflows/local/phylogeny/main'
include { SNPEFF               } from '../subworkflows/local/snpeff/main'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { paramsSummaryMap     } from 'plugin/nf-schema'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/fastqc/main'
include { QC_REPORTSHEET              } from '../modules/local/qc_reportsheet/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { GATK4_HAPLOTYPECALLER       } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { SEQKIT_REPLACE              } from '../modules/nf-core/seqkit/replace/main'
include { SNPDISTS                    } from '../modules/nf-core/snpdists/main'
include { GATK4_COMBINEGVCFS          } from '../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_INDEXFEATUREFILE      } from '../modules/nf-core/gatk4/indexfeaturefile/main'
include { QC_PARSER                   } from '../modules/local/qc_parser/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary

workflow MYCOSNP {

    take:
    samplesheet // New samplesheet combines ingestion for fastq reads, sra accessions, and vcf files for phylogeny

    main:
    if ( !params.fasta && !params.ref_dir && !params.ref_masked_fasta && !params.ref_fai && !params.ref_bwa && !params.ref_dict) {
        log.error "Genome fasta or index files not specified with e.g. '--fasta genome.fa', '--ref_dir genome/' or via a detectable config file."
        System.exit(1)
    }

    def ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        params.input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_all_reads    = INPUT_CHECK.out.ch_fastq_reads  // channel: [ val(meta), [ reads ] ]
    ch_sra_list     = INPUT_CHECK.out.ch_sra_list     // channel: [ val(meta), sra_id    ]
    ch_vcf_files    = INPUT_CHECK.out.ch_vcf_files    // channel: [ val(meta), vcf       ]
    ch_vcfidx_files = INPUT_CHECK.out.ch_vcfidx_files // channel: [ val(meta), tbi       ]

    //
    // SUBWORKFLOW: Fetch FASTQ reads from input SRA accession IDs, mix with fastq reads from samplesheet
    //
    SRA_FASTQ_SRATOOLS (
        ch_sra_list
    )
    ch_versions  = ch_versions.mix(SRA_FASTQ_SRATOOLS.out.versions)

    ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)

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

    if(params.ref_dir != null) {
        def meta_ref = [ id: 'reference' ]
        ch_ref_fasta = Channel.fromPath(params.ref_dir + "/masked/*.fa*", checkIfExists:true)
                             .first()
                             .map { p -> [ meta_ref, p ] }
        ch_ref_fai   = Channel.fromPath(params.ref_dir + "/fai/*.fai", checkIfExists:true)
                             .first()
                             .map { p -> [ meta_ref, p ] }
        ch_ref_bwa   = Channel.fromPath(params.ref_dir + "/bwa/bwa", checkIfExists:true, type: 'dir')
                             .first()
                             .map { p -> [ meta_ref, p ] }
        ch_ref_dict  = Channel.fromPath(params.ref_dir + "/dict/*.dict", checkIfExists:true)
                             .first()
                             .map { p -> [ meta_ref, p ] }

    } else if ( params.ref_masked_fasta && params.ref_fai && params.ref_bwa && params.ref_dict ) {
        def meta_ref = [ id: 'reference' ]
        if(params.ref_masked_fasta != null) {
            ch_ref_fasta  = Channel.fromPath(params.ref_masked_fasta, checkIfExists:true).first().map { p -> [ meta_ref, p ] }
        }

        if(params.ref_fai != null) {
            ch_ref_fai    = Channel.fromPath(params.ref_fai, checkIfExists:true).first().map { p -> [ meta_ref, p ] }
        }

        if(params.ref_bwa != null) {
            ch_ref_bwa    = Channel.fromPath(params.ref_bwa, checkIfExists:true, type: 'dir').first().map { p -> [ meta_ref, p ] }
        }

        if(params.ref_dict != null) {
            ch_ref_dict   = Channel.fromPath(params.ref_dict, checkIfExists:true).first().map { p -> [ meta_ref, p ] }
        }
    } else if (params.fasta) {
        ch_fasta = file(params.fasta)
        BWA_REFERENCE (
            ch_fasta,
            params.mask
        )
        ch_versions = ch_versions.mix(BWA_REFERENCE.out.versions)

        // Use tuple channels directly from reference subworkflow
        ch_ref_fasta = BWA_REFERENCE.out.masked_fasta
        ch_ref_fai   = BWA_REFERENCE.out.samtools_index
        ch_ref_bwa   = BWA_REFERENCE.out.bwa_index
        ch_ref_dict  = BWA_REFERENCE.out.dict
    } else {
        exit 1, 'Input reference fasta or index files not specified!'
    }

    // Derive path-only channels for modules that expect plain file inputs
    ref_fasta_only = ch_ref_fasta.map{ meta1, fa1 -> fa1 }.first()
    ref_fai_only   = ch_ref_fai.map{ meta1, fai -> fai }.first()
    ref_dict_only  = ch_ref_dict.map{ meta1, dict -> dict }.first()

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
    // Pass reference channels separately to BWA_PREPROCESS
    BWA_PREPROCESS (
        ch_ref_fasta,
        ch_ref_fai,
        ch_ref_bwa,
        ch_all_reads,
        params.coverage
    )
    ch_versions = ch_versions.mix(BWA_PREPROCESS.out.versions)

    // MODULE: QC_REPORTSHEET
    ch_qcreportsheet = BWA_PREPROCESS.out.qc_lines.collect()
    QC_REPORTSHEET (
        ch_qcreportsheet
    )

    // Conditionally run QC_PARSER if param.amdp is true
    if (params.amdp) {
        QC_PARSER (
            QC_REPORTSHEET.out.qc_reportsheet
        )
        ch_versions  = ch_versions.mix(QC_PARSER.out.versions)
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

    // Pad alignment tuple with placeholders for intervals and dragstr_model expected by the module's first input port
    def hc_align = BWA_PREPROCESS.out.alignment_combined.map { meta, bam, bai -> [ meta, bam, bai, [], [] ] }
    GATK4_HAPLOTYPECALLER (
        hc_align,
        ch_ref_fasta,
        ch_ref_fai,
        ch_ref_dict,
        [[],[]],  // dbsnp (optional)
        [[],[]]   // dbsnp_tbi (optional)
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    if(! params.skip_combined_analysis) {

        def combined_meta_id = [ id:'combined', single_end:false ]

        // vcf files
        ch_vcf_hc = GATK4_HAPLOTYPECALLER.out.vcf.map{ meta, vcf -> vcf }
        ch_vcf_raw = ch_vcf_files.map { meta, vcf -> vcf }
        ch_vcf_all = ch_vcf_hc.concat(ch_vcf_raw).toSortedList().map { vcfs -> [ combined_meta_id, vcfs ] }
        // index files
        ch_tbi_hc = GATK4_HAPLOTYPECALLER.out.tbi.map{ meta, tbi -> tbi }
        ch_tbi_raw = ch_vcfidx_files.map { meta, tbi -> tbi }
        ch_tbi_all = ch_tbi_hc.concat(ch_tbi_raw).toSortedList().map { tbi -> [ combined_meta_id, tbi ] }

        GATK4_COMBINEGVCFS (
            ch_vcf_all.join(ch_tbi_all),
            ref_fasta_only,
            ref_fai_only,
            ref_dict_only
        )

        ch_versions = ch_versions.mix(GATK4_COMBINEGVCFS.out.versions)

        GATK4_INDEXFEATUREFILE (
            GATK4_COMBINEGVCFS.out.combined_gvcf
        )

        ch_versions = ch_versions.mix(GATK4_INDEXFEATUREFILE.out.versions)

        ch_gatk_variants = GATK4_COMBINEGVCFS.out.combined_gvcf.join( GATK4_INDEXFEATUREFILE.out.index )

        GATK_VARIANTS (
            ref_fasta_only,
            ref_fai_only,
            ch_ref_bwa.map { m, b -> b }.first(),
            ref_dict_only,
            ch_gatk_variants,
            params.max_amb_samples,
            params.max_perc_amb_samples,
            params.min_depth

        )
        ch_versions = ch_versions.mix(GATK_VARIANTS.out.versions)

        if(params.snpeff != false) {
            SNPEFF (
                GATK_VARIANTS.out.filtered_vcf,
                params.species,
                params.snpeffcache,
                params.ref_name
            )
            ch_versions = ch_versions.mix(SNPEFF.out.versions)
        }

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
        SEQKIT_REPLACE (
            GATK_VARIANTS.out.snps_fasta
        ) // Swap * for -
        ch_versions = ch_versions.mix(SEQKIT_REPLACE.out.versions)

        SNPDISTS (
            SEQKIT_REPLACE.out.fastx
        )
        ch_versions = ch_versions.mix(SNPDISTS.out.versions)

        if(!params.skip_phylogeny) {
            CREATE_PHYLOGENY (
                SEQKIT_REPLACE.out.fastx.map{meta, fas->[fas]},
                '',
                SNPDISTS.out.tsv,
                params.rapidnj,
                params.fasttree,
                params.iqtree,
                params.raxmlng,
                params.amdp,
                params.metadata_csv,
                params.geolocation_csv
            )
            ch_versions = ch_versions.mix(CREATE_PHYLOGENY.out.versions)
        }

    }

    //
    // MODULE: Run Pre-FastQC
    //
    FASTQC_RAW (
        ch_all_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'mycosnp_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    workflow_summary    = paramsSummaryMultiqc(summary_params)
    ch_workflow_summary = Channel.from(workflow_summary)

    def ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.post_qc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.stats.map{meta, stats -> [stats]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.flagstat.map{meta, stats -> [stats]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.idxstats.map{meta, stats -> [stats]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.qualimap.map{meta, stats -> [stats]}.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        [],
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report
    versions = ch_versions

}
/*
========================================================================================
    THE END
========================================================================================
*/
