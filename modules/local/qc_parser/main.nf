process QC_PARSER {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    path qc_reportsheet

    output:
    path("combined_QC_results.csv"), emit: qc_fail_pass

    script:
    def args = task.ext.args ?: ''
    def qc_thresholds = params.qc_thresholds

    """
    qc_parser.py ${qc_reportsheet} -qc_thresholds ${qc_thresholds} \\
        ${args}
    """

    stub:
    """
    touch combined_QC_results.csv
    """
}
