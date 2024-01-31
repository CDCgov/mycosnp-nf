process SOURMASH_SKETCH {
    tag "$meta.id"
    label 'process_single'
    //copied from nf-core modules and modified

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.4--hdfd78af_0':
        'biocontainers/sourmash:4.8.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("*.sig"), emit: signatures
    path("*.sig")                 , emit: sig_paths
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // required defaults for the tool to run, but can be overridden
    def args = task.ext.args ?: "dna --param-string 'k=31'"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash sketch \\
        $args \\
        --output '${prefix}.sig' \\
        $sequence

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}