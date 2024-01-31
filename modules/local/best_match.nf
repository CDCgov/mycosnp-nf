process BEST_MATCH {
    tag "$meta.id"
    label 'process_single'
    //takes modified matrix and removes non-sample comparisons, sorts by best match
    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' : 
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(mod_matrix)

    output:
    tuple val(meta), path('*.csv')      , emit: indv_bestmatch
    path('*.csv')                       , emit: match_files
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/mycosnp/bin/
    def prefix         = task.ext.prefix ?: "${meta.id}"
    """
    best_match.py \\
     $mod_matrix \\
     $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
