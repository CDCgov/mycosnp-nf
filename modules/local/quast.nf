process QUAST {
    label 'process_medium'

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    tuple val(meta), path(ref), path(reads), path(assembly) 

    output:
    path "${prefix}"    , emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "$meta.id"
    """
    quast.py \\
        -r ${ref} \\
        --output-dir ${prefix} \\
        --threads $task.cpus \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        $args \\
        $assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
