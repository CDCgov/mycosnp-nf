process INPUT_PROC {

    conda (params.enable_conda ? "conda-forge::python=3.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
    containerOptions { workflow.containerEngine == "docker" ? '--tmpfs /tmp:rw,noexec,nosuid,size=1g --tmpfs /var/tmp:rw,noexec,nosuid,size=1g' : '' }

    input:
    path(fasta)

    output:
    tuple val(meta), path("reference.fasta"), path("reference.copy.fasta"), emit: ref_fasta
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta = ["id": 'reference']
    def is_compressed = false
    if(fasta.getName().endsWith(".gz"))
    {
        is_compressed = true
    }
    """
    echo ${meta.id}
    if [[ ${is_compressed} == "true" ]]; then
        gunzip -c $fasta > reference.fasta
    else
        cp ${fasta} reference.fasta
    fi
    cp reference.fasta reference.copy.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    meta = ["id": 'reference']
    """
    touch reference.fasta
    touch reference.copy.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
