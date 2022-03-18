//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { LANE_MERGE        } from '../../modules/local/lane_merge'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    Channel.fromPath(samplesheet)
           .splitCsv( header:false, sep:',', skip:1 )
           .map { row -> stage_fastq(row) }
           .set{ precheck_reads }

    LANE_MERGE(precheck_reads)

    emit:
    reads    =   LANE_MERGE.out.reads       // channel: [ val(meta), [ reads ] ]
    versions =   ch_versions                // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2, fastq_?... ] ]
def stage_fastq(ArrayList row) {
    //print row
    def meta        = [:]
    meta.id         = row[0]
    meta.single_end = false
    def array       = []
    def filesarray  = []

    for(int i = 1; i < row.size(); i++)
    {
        if(row[i] == "")
        {
            // skip this row
        } else if (!file(row[i]).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read $i FastQ file does not exist!\n${row[i]}"
        } else
        {
            filesarray.add(file(row[i]))
        }
    }

    if(filesarray.size() == 1)
    {
        meta.single_end = true
    } else if( (filesarray.size() % 2) != 0)
    {
        exit 1, "ERROR: Please check input samplesheet -> Number of samples is not an even number or 1.\n$row"
    }
    
    array = [ meta, filesarray]
    return array
}
