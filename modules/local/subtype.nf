process SUBTYPE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mash=2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash%3A2.3--he348c14_1' : 
        'quay.io/staphb/mash:2.3' }"
        
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
    species=\$(cat ${gambit_results} | grep -v "predicted.name" | tr ' ' '_' | tr ',' '\t' | cut -f 2)
    
    # run subtyper
    subtyper.sh ${prefix} \${species} ${subtype_db} ${seq} || true

    # check if file was created, if not then create empty file
    if [ ! -f "${prefix}_subtype.csv" ]
    then
        touch ${prefix}_subtype.csv
    fi
    """
}
