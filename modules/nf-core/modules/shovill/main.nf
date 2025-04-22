process SHOVILL {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::shovill=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shovill:1.1.0--0' :
        'biocontainers/shovill:1.1.0--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("contigs.fa")                         , emit: contigs
    tuple val(meta), path("shovill.log")                        , emit: log
    tuple val(meta), path("{skesa,spades,megahit,velvet}.fasta"), emit: raw_contigs
    tuple val(meta), path("contigs.{fastg,gfa,LastGraph}")      , optional:true, emit: gfa
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gsize = params.genome_size ? "--gsize ${params.genome_size}" : ''
    def args = task.ext.args ?: ''
    def memory = task.memory.toGiga()
    """
    mkdir -p tmp
    shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        $args \\
        $gsize \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./ \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(echo \$(shovill --version 2>&1) | sed 's/^.*shovill //')
    END_VERSIONS
    """
}
