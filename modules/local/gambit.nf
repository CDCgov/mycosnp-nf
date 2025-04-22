process GAMBIT_QUERY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gambit=1.0.0" : null)
    container 'jdj0303/gambit:1.0.0'
    
    input:
    tuple val(meta), path(assembly)
    path db_file
    path h5_file

    output:
    tuple val(meta), path("*_gambit.csv"), emit: taxa
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gambit \\
        ${args} \\
        -d ./ \\
        query ${assembly} > ${prefix}_gambit.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gambit: \$(gambit --version | tr -d "[a-z, ]")
    END_VERSIONS
    """
}
