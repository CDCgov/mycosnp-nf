process PRE_MYCOSNP_COMB_SUMMARY {
    label 'process_low'
    
    container 'ubuntu:jammy'

    input:
    path line_summary

    output:
    path "pre-mycosnp-summary.csv"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # combine all line summaries into single report
    echo "Sample,Predicted_Rank,Predicted_Taxon,Subtype_Clade,Subtype_Reference_Accession,Subtype_ANI,Closest_GAMBIT_Entry_Description,Closest_GAMBIT_Entry_Distance,Trimmed_Reads,Avg_Read_Quality,Sample_Assembly_Length,Sample_Assembly_GC,Reference_Genome_Length,Avg_Depth_Coverage,Reference_GC" > pre-mycosnp-summary.csv
    cat *_linesummary.csv >> pre-mycosnp-summary.csv 
    """
}
