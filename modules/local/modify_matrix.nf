process MODIFY_MATRIX {
    tag "$meta.id"
    label 'process_single'

    //modifies an ANI matrix to a simple network file
    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' : 
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(ani_matrix)

    output:
    tuple val(meta), path('*.tsv')       , emit: tsv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/mycosnp/bin/
    def prefix         = task.ext.prefix ?: "${meta.id}"
    """
    modmatrix.py \\
        $ani_matrix \\
        ${prefix}.ani_comp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
