process PRE_MYCOSNP_COMB_SUMMARY {
    label 'process_low'

    conda "conda-forge::coreutils=9.5"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
    'quay.io/biocontainers/coreutils:9.5'}"

    input:
    path line_summary

    output:
    path "pre-mycosnp-summary.csv"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # combine all line summaries into single report
    echo "Sample,PM_Predicted_Rank,PM_Predicted_Taxon,PM_Subtype_Closest_Match,PM_Subtype_ANI,PM_Closest_GAMBIT_Entry_Description,PM_Closest_GAMBIT_Entry_Distance,PM_Trimmed_Reads,PM_Avg_Read_Quality,PM_Sample_Assembly_Length,PM_Sample_Assembly_GC,PM_Reference_Genome_Length,PM_Avg_Depth_Coverage,PM_Reference_GC" > pre-mycosnp-summary.csv
    cat *_linesummary.csv >> pre-mycosnp-summary.csv
    """
}
