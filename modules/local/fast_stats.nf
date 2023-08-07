process FAST_STATS {
    label 'process_low'
    
    container 'ubuntu:jammy'   

    input:
    tuple val(meta), path(assembly), path(ref), path(faqcs) 

    output:
    tuple val(meta), path("*_fast-stats.csv")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "$meta.id"
    """
    # run fast-stats
    fast-stats.sh \\
        ${prefix} \\
        ${assembly} \\
        ${ref} \\
        *.stats.txt \\
        *.for_qual_histogram.txt \\
        > ${prefix}_fast-stats.csv
    
    """
}
