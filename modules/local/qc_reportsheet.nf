process QC_REPORTSHEET {
    label 'process_low'

    input:
    path(qc_lines)

    output:
    path("qc_reportsheet.csv"), emit: qc_reportsheet

    script:

    """
    printf \"Sample Name\\t# Reads Before Trimming\\tGC Before Trimming\\tAverage Phred Before Trimming\\tCoverage Before Trimming\\t# Reads After Trimming\\t# Paired Reads After Trimming\\t# Unpaired Reads After Trimming\\tGC After Trimming\\tAverage Phred After Trimming\\tCoverage After Trimming\\n\" > qc_reportsheet.csv
    cat ${qc_lines} >> qc_reportsheet.csv
    """

}

