process PICARD_ADDORREPLACEREADGROUPS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args        ?: ''
    def prefix = task.ext.prefix    ?: "${meta.id}"
    def ID = task.ext.id            ?: "id"
    def LIBRARY= task.ext.library   ?: "library"
    def PLATFORM= task.ext.platform ?: "illumina"
    def BARCODE= task.ext.barcode   ?: "barcode"
    def SAMPLE= task.ext.sample     ?: "sample"
    def INDEX= task.ext.index       ?: "index"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard AddOrReplaceReadGroups] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        AddOrReplaceReadGroups \\
        -Xmx${avail_mem}g \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.bam \\
        -ID ${ID} \\
        -LB ${LIBRARY} \\
        -PL ${PLATFORM} \\
        -PU ${BARCODE} \\
        -SM ${SAMPLE} \\
        -CREATE_INDEX true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard AddOrReplaceReadGroups --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix    ?: "${meta.id}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard AddOrReplaceReadGroups --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
