process SUBTYPE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.14--hdfd78af_0' :
        'biocontainers/sourmash:4.8.14--hdfd78af_0' }"

    input:
    tuple val(meta), path(gambit_results), path(seq)
    path(subtype_db)

    output:
    tuple val(meta), path("*_subtype.csv"), emit: subtype
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # determine species call from Gambit
    taxon=\$(cat ${gambit_results} | grep -v "predicted.name" | tr ' ' '_' | cut -f 2 -d ',')

    # run subtyper only if species is not blank
    if [ -n "\$taxon" ]; then
        subtyper.py ${prefix} \${taxon} ${subtype_db} ${seq} ${prefix}_subtype.csv
    fi

    # check if file was created, if not then create empty file
    if [ ! -f "${prefix}_subtype.csv" ]
    then
        touch ${prefix}_subtype.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | tr -d "[a-z, ]")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_subtype.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | tr -d "[a-z, ]")
    END_VERSIONS
    """
}
