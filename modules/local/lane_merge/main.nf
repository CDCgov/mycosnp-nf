process LANE_MERGE {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("combined/*.fastq.gz"), emit: reads

    script:
    def reads_list = reads.collect { "'${it.name}'" }.join(' ')
    """
    mkdir -p combined

    READS=(${reads_list})
    NUM_READS=\${#READS[@]}

    if [[ "\$NUM_READS" -eq 0 ]]; then
        echo "[COMBINE_READS] No fastqs provided for sample: ${meta.id}" >&2
        exit 1
    fi

    echo "[COMBINE_READS] Sample: ${meta.id} | Files (\$NUM_READS): \${READS[*]}" >&2

    if [[ "\$NUM_READS" -eq 1 ]]; then
        f="\${READS[0]}"
        if [[ "\$f" == *.gz ]]; then
            zcat "\$f" | gzip -c > "combined/${meta.id}.fastq.gz"
        else
            gzip -c "\$f" > "combined/${meta.id}.fastq.gz"
        fi
    else
        # Split by index cardinality into R1 (0,2,4,...) and R2 (1,3,5,...) assuming given order corresponds to read pair mates
        R1=()
        R2=()
        for i in "\${!READS[@]}"; do
            if (( i % 2 == 0 )); then
                R1+=( "\${READS[\$i]}" )
            else
                R2+=( "\${READS[\$i]}" )
            fi
        done

        # Combine R1 (assumes files follow expected file cardinality)
        {
            for f in "\${R1[@]}"; do
                if [[ "\$f" == *.gz ]]; then
                    zcat "\$f"
                else
                    cat "\$f"
                fi
            done
        } | gzip -c > "combined/${meta.id}_R1.fastq.gz"

        # Combine R2 files (again assumes file cardinality structuring)
        {
            for f in "\${R2[@]}"; do
                if [[ "\$f" == *.gz ]]; then
                    zcat "\$f"
                else
                    cat "\$f"
                fi
            done
        } | gzip -c > "combined/${meta.id}_R2.fastq.gz"
    fi
    """
}
