process COORDSTOBED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mummer=3.23" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mummer:3.23--pl5262h1b792b2_12' :
        'quay.io/biocontainers/mummer:3.23--pl5262h1b792b2_12' }"

    input:
    tuple val(meta), path(delta)

    output:
    tuple val(meta), path("masked_ref.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "masked_ref"
    """
    show-coords ${args} $delta > masked_ref_BEFORE_ORDER.bed

    awk '{if (\$1 != \$3 && \$2 != \$4) print \$0}' masked_ref_BEFORE_ORDER.bed > masked_ref_BEFORE_ORDER2.bed
    
    awk '{print \$8\"\\t\"\$1\"\\t\"\$2}' masked_ref_BEFORE_ORDER2.bed > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mummer: \$(nucmer --version 2>&1 | grep "version" | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch masked_ref.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mummer: \$(nucmer --version 2>&1 | grep "version" | sed 's/^.*version //')

    END_VERSIONS
    """
}
