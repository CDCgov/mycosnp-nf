process SUBTYPE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sourmash-minimal=4.8.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.14--hdfd78af_0' :
        'quay.io/biocontainers/sourmash:4.8.14--hdfd78af_0' }"
        
    input:
    tuple val(meta), path(gambit_results), path(seq)
    path subtype_db

    output:
    tuple val(meta), path("*_subtype.csv"), emit: subtype

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
    """
}
