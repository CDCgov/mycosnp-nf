process REPORT_SOURMASH {
    //tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' : 
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    path(match_files)

    output:
    path('*.csv')      , emit: final_report
    
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/mycosnp/bin/
    //def prefix         = task.ext.prefix ?: "${meta.id}"
    // _mqc.csv
    """
    sourmash_report.py \\
    $match_files
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
