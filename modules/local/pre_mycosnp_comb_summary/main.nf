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
    echo "Sample,Predicted_Rank,Predicted_Taxon,Subtype_Closest_Match,Subtype_ANI,Closest_GAMBIT_Entry_Description,Closest_GAMBIT_Entry_Distance,Trimmed_Reads,Avg_Read_Quality,Sample_Assembly_Length,Sample_Assembly_GC,Reference_Genome_Length,Avg_Depth_Coverage,Reference_GC" > pre-mycosnp-summary.csv
    cat *_linesummary.csv >> pre-mycosnp-summary.csv
    """
}
