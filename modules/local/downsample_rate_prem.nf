process DOWNSAMPLE_RATE_PREM {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
	
    input:
        tuple val(meta), path(reads)
        path(reference_fasta)
        val(coverage)

    output:
        tuple val(meta), env(SAMPLE_RATE),  env(SAMPLED_NUM_READS),   emit: downsampled_rate
	env(SAMPLED_NUM_READS)                                    ,   emit: number_to_sample

	script:
	"""
	REFERENCE_LEN=\$(awk '!/^>/ {len+=length(\$0)} END {print len}' < ${reference_fasta})
	READS_LEN=\$(zcat ${reads} |awk 'NR%4==2 {len +=length(\$0)} END {print len}')
	
	SAMPLE_RATE=\$(echo "${coverage} \${READS_LEN} \${REFERENCE_LEN}" | awk '{x=\$1/(\$2/\$3); x=(1<x?1:x)} END {print x}')
	
	# Calculate number of reads
	NUM_READS=\$(zcat ${reads[0]}|awk 'END {print NR/4}')	
	SAMPLED_NUM_READS=\$(echo "\${NUM_READS} \${SAMPLE_RATE}" | awk '{x=\$1*\$2} END {printf "%.0f", x}')
	"""
}
