/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMycosnp.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
if (params.skip_samples_file) { // check for skip_samples_file
    checkPathParamList.add(params.skip_samples_file)
}
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
sra_list = []
sra_ids = [:]
ch_input = null;
if (params.input) 
{ 
    ch_input = file(params.input) 
}
if(params.add_sra_file)
{
    sra_file = file(params.add_sra_file, checkIfExists: true)
    allLines  = sra_file.readLines()
    for( line : allLines ) 
    {
        row = line.split(',')
        if(row.size() > 1)
        {
            println "Add SRA ${row[1]} => ${row[0]}"
            sra_list.add(row[1])
            sra_ids[row[1]] = row[0]
        } else
        {
            if(row[0] != "")
            {
                println " ${row[0]} => ${row[0]}"
                sra_list.add(row[0])
                sra_ids[row[0]] = row[0]
            }
        }
    }
}

vcf_file_list = []
vcfidx_file_list = []
if(params.add_vcf_file)
{
    vcf_file = file(params.add_vcf_file, checkIfExists: true)
    allLines  = vcf_file.readLines()
    for( line : allLines ) 
    {
        if(line != "")
        {
            println " Add VCF => $line"
            t_vcf = file(line)
            t_idx = file(line + ".tbi")
            vcf_file_list.add(t_vcf)
            vcfidx_file_list.add(t_idx)
        }
    }
}

if(! ( params.input || params.add_sra_file || params.add_vcf_file ) ) { exit 1, 'Input samplesheet, sra file, or vcf file not specified!' }


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
include { SRA_FASTQ_SRATOOLS } from '../subworkflows/local/sra_fastq_sratools'
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { BWA_PREPROCESS     } from '../subworkflows/local/bwa-pre-process'
include { BWA_REFERENCE      } from '../subworkflows/local/bwa-reference'
include { GATK_VARIANTS      } from '../subworkflows/local/gatk-variants'
include { CREATE_PHYLOGENY   } from '../subworkflows/local/phylogeny'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/modules/fastqc/main'
include { QC_REPORTSHEET              } from '../modules/local/qc_reportsheet.nf'
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
    ch_all_reads = Channel.empty()
    ch_sra_reads = Channel.empty()
    ch_sra_list  = Channel.empty()
    if(params.add_sra_file)
    {   
        ch_sra_list = Channel.fromList(sra_list)
                             .map{valid -> [ ['id':sra_ids[valid],single_end:false], valid ]}
        SRA_FASTQ_SRATOOLS(ch_sra_list)
        ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)
    }
    
    if(params.input)
    {
        INPUT_CHECK (
            ch_input
        )
        ch_all_reads = ch_all_reads.mix(INPUT_CHECK.out.reads)
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    }


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

    BWA_PREPROCESS( [fas_file, fai_file, bai_file ], ch_all_reads)
    ch_versions = ch_versions.mix(BWA_PREPROCESS.out.versions)

    // MODULE: QC_REPORTSHEET
    ch_qcreportsheet = BWA_PREPROCESS.out.qc_lines.collect()
    QC_REPORTSHEET (
        ch_qcreportsheet
    )

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

    GATK4_HAPLOTYPECALLER(  BWA_PREPROCESS.out.alignment_combined.map{meta, bam, bai            -> [ meta, bam, bai, [] ] },
                            fas_file,
                            fai_file,
                            dict_file,
                            [],
                            []
                            )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)


    ch_vcf_files = Channel.empty()
    if(! params.skip_combined_analysis)
    {

        ch_vcf_files    = Channel.fromList(vcf_file_list)
        ch_vcfidx_files = Channel.fromList(vcfidx_file_list)

        ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{meta, vcf ->[ vcf ]  }.collect()
        ch_vcf = ch_vcf.mix(ch_vcf_files)
        ch_vcf_idx = GATK4_HAPLOTYPECALLER.out.tbi.map{meta, idx ->[ idx ]  }.collect()
        ch_vcf_idx = ch_vcf_idx.mix(ch_vcfidx_files)
        
        ch_vcfs = ch_vcf.mix(ch_vcf_idx).collect()

        GATK4_LOCALCOMBINEGVCFS(
                                    [id:'combined', single_end:false],
                                    ch_vcfs,
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
        if(! params.skip_phylogeny) {
            CREATE_PHYLOGENY(SEQKIT_REPLACE.out.fastx.map{meta, fas->[fas]}, '')
        }
    }

     CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: Run Pre-FastQC 
    //
    FASTQC_RAW (
        ch_all_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BWA_PREPROCESS.out.post_qc.collect{it[1]}.ifEmpty([]))
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
