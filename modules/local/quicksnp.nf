process QUICKSNP {
    tag "$meta.id"
    label 'process_low'

    // conda (params.enable_conda ? "bioconda::quicksnp=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/staphb/quicksnp:1.0.1' :
        'quay.io/staphb/quicksnp:1.0.1' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.nwk"), emit: quicksnp_tree
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    QuickSNP.py \\
        --dm ${tsv} \\
        --outtree quicksnp_tree.nwk
    """
}