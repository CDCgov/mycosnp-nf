process QC_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::scipy=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scipy%3A1.1.0' :
        'quay.io/biocontainers/scipy:1.1.0' }"
        
    input:
    tuple val(meta), path(txt)
    //path reference

    output:
    tuple val(meta), path("output.txt"), emit: qc_stuff

    """
    python $projectDir/bin/qc_report_stats.py > output.txt

    """
}