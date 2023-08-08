process COMBINE_SUMMARY {
    label 'process_low'
    
    container 'ubuntu:jammy'

    input:
    path line_summary

    output:
    tuple path("pre-mycosnp-summary.csv")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # combine all line summaries into single report
    echo "Sample,Trimmed_Reads,Avg_Read_Quality,Avg_Depth_Coverage,Sample_Genome_Length,Reference_Genome_Length,Sample_GC,Reference_GC,Predicted_Taxa,Predicted_Subtype,Subtype_ANI,Reference_Accession" > pre-mycosnp-summary.csv
    cat *_linesummary.csv >> pre-mycosnp-summary.csv 
    """
}
