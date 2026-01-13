//
// Check input samplesheet and get read channels
//

include { LANE_MERGE        } from '../../../modules/local/lane_merge/main'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    // "Normalize" each row to triage out into channels
    // for each file type downstream, keep only non-empty strings.
    rows_ch = samplesheet
        .splitCsv(header:true, sep: ',')
        .map { row ->
            def norm = { col -> col == null ? null : col.toString().trim() }
            def record = [
                sample : norm(row.sample),
                fastq1 : norm(row.fastq_1),
                fastq2 : norm(row.fastq_2),
                sra    : norm(row.sra),
                vcf    : norm(row.vcf)
            ]
            if( !record.sample ) {
                exit 1, "ERROR: Sample name is required for all samples. Please check your input samplesheet."
            }
            return record
        }

    // Triage checkpoint 1: Split out a channel of fastq reads --> results in tuple(meta, [reads])
    // Includes a checkpoint to ensure that all fastq files exist.
    fastq_pairs_flat = rows_ch
        .flatMap { row ->
            [row.fastq1, row.fastq2, row.fastq3, row.fastq4]
            .findAll { it }
            .collect { fq ->
                def f = file(fq)
                if( !f.exists() )
                    exit 1, "ERROR: Please check input samplesheet -> FastQ file does not exist!\n${fq} (sample: ${row.sample})"
                tuple(row.sample, f)
            }
        }

    // break up read sets into multi-lane reads and single-/paired-end reads
    precheck_reads = fastq_pairs_flat
        .groupTuple(by: 0)
        .map { sample, fqs_list ->
            List fqs = (fqs_list as List)
            def uniq_reads = fqs.unique(false).sort { it.getName() }                // Sort should preserve the ordering of fastq files here for LANE_MERGE downstream
            def meta = [
                id: sample,
                single_end: uniq_reads.size() == 1                                  // I don't think we're technically using this, but kept it for consistency
            ]
            tuple(meta, uniq_reads)
        }
        .branch {
            multi_lane: it[1].size() > 2 // filter samples with > 2 fastq files
            single_pair: it[1].size() <= 2 // filter single-/paired-end files
        }

    //only process samples with mutiple lanes
    LANE_MERGE (
        precheck_reads.multi_lane
    )
    ch_versions = ch_versions.mix(LANE_MERGE.out.versions)

    // Triage checkpoint 2: Split out a channel of SRA accessions to fetch --> results in tuple(meta, sra_id)
    ch_sra_list = rows_ch
        .map { row ->
            if( row.sra ) {
                def meta = [
                    id        : row.sample,
                    single_end: false               // Sanity check: Assume paired by default for SRA
                ]
                tuple(meta, row.sra)
            } else {
                null
            }
        }
        .filter { it != null }

    // Triage checkpoint 3: Split out a channel of VCF files for phylogeny downwind --> results in tuple(meta, vcf_file)
    // Includes logic to ensure that VCF files exist and that it has a corresponding .tbi index file.
    vcf_triples = rows_ch
        .map { row ->
            // Ignore if VCF column not filled out
            if( !row.vcf ) {
                return null
            }

            // If VCF column filled out, check that the file exists
            def vcf = file(row.vcf)
            if( !vcf.exists() ) {
                exit 1, "ERROR: VCF file not found: ${row.vcf} (sample: ${row.sample})"
            }

            // If checks pass so far, also check that the tabix index "<vcf>.tbi" also exists
            // This should handle both local files as well as S3
            def tbi = vcf.getParent().resolve(vcf.getName() + '.tbi')
            if( !tbi.exists() ) {
                exit 1, "ERROR: Tabix index (.tbi) not found for VCF: ${vcf} (sample: ${row.sample})\n" +
                    "Expected index path: ${tbi}\n" +
                    "Tip: bgzip + tabix (e.g., 'bgzip -c sample.vcf > sample.vcf.gz' then 'tabix -p vcf sample.vcf.gz')."
            }

            // Okay, all file checks passed.  Build the tuple.
            def meta = [
                id        : row.sample,
                single_end: null                // not applicable for VCFs
            ]
            tuple(meta, vcf, tbi)
        }
        .filter { it != null }

        // Split the 3-element tuple into two channels, one for vcfs and one for tbi indices
        ch_vcf_idx    = vcf_triples.map { meta, vcf, tbi -> tuple(meta, tbi) }
        ch_vcf_files  = vcf_triples.map { meta, vcf, tbi -> tuple(meta, vcf) }


    emit:
    ch_fastq_reads  =   precheck_reads.single_pair.mix( LANE_MERGE.out.reads )      // channel: [ val(meta), [ reads ] ]
    ch_sra_list     =   ch_sra_list                                                 // channel: [ val(meta), sra_id    ]
    ch_vcf_files    =   ch_vcf_files                                                // channel: [ val(meta), vcf       ]
    ch_vcfidx_files =   ch_vcf_idx                                                  // channel: [ val(meta), tbi       ]
    versions        =   ch_versions                                                 // channel: [ versions.yml         ]
}
