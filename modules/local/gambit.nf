process GAMBIT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gambit=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gambit%3A1.0.0--py39hbf8eff0_0' : 
        'quay.io/staphb/gambit:1.0.0' }"
        
    input:
    tuple val(meta), path(assembly)
    path db_dir

    output:
    path("*_gambit.txt"), emit: taxa

    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gambit \\
        ${args} \\
        -d ${db_dir} \\
        query ${assembly} > ${prefix}_gambit.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gambit: \$(gambit --version | tr -d "[a-z, ]")
    END_VERSIONS
    """
}
