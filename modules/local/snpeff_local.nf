process SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::snpeff=4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:4.3.1t--hdfd78af_5' :
        'quay.io/biocontainers/snpeff:4.3.1t--hdfd78af_5' }"

    input:
    tuple val(meta), path(vcf)
    val species
    path(config)
    path(datadir)


    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    path "*.csv"                      , emit: report
    path "*.html"                     , emit: summary_html
    path "*.genes.txt"                , emit: genes_txt
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -r $datadir data

    cp $config snpEff.config

    snpEff \\
        -Xmx${avail_mem}g \\
        -config snpEff.config \\
        -v $species \\
        -csvStats ${prefix}.csv \\
        $vcf \\
        > ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}