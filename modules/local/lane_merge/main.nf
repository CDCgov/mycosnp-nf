process LANE_MERGE {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads

    script:
    def reads_list = reads.collect { "'${it.name}'" }.join(' ')
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    READS=(${reads_list})
    NUM_READS=\${#READS[@]}

    # Handles samples with mutiple lanes (sequences) and concatenates
    # Split by index cardinality into R1 (0,2,4,...) and R2 (1,3,5,...)

    R1=()
    R2=()

    for i in "\${!READS[@]}"; do
        if (( i % 2 == 0 )); then
            R1+=( "\${READS[\$i]}" )
        else
            R2+=( "\${READS[\$i]}" )
        fi
    done

    # Process R1 and R2 in parallel using pigz
    {
        for f in "\${R1[@]}"; do
            if [[ "\$f" == *.gz ]]; then
                zcat "\$f"
            else
                cat "\$f"
            fi
        done
    } | pigz -p ${task.cpus} > "${prefix}_R1.fastq.gz" &

    {
        for f in "\${R2[@]}"; do
            if [[ "\$f" == *.gz ]]; then
                zcat "\$f"
            else
                cat "\$f"
            fi
        done
    } | pigz -p ${task.cpus} > "${prefix}_R2.fastq.gz" &

    wait
    """

    stub:
    def reads_list = reads.collect { "'${it.name}'" }.join(' ')
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    READS=(${reads_list})
    NUM_READS=\${#READS[@]}

    if [[ "\$NUM_READS" -eq 1 ]]; then
        # Single file case - create empty gzipped file
        echo -n | pigz -p ${task.cpus} > "${prefix}.fastq.gz"
    else
        # Multiple files case - assume paired-end
        echo -n | pigz -p ${task.cpus} > "${prefix}_R1.fastq.gz"
        echo -n | pigz -p ${task.cpus} > "${prefix}_R2.fastq.gz"
    fi
    """
}
