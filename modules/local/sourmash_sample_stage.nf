process SOURMASH_SAMPLE_STAGE {
    tag "$meta.id"
    label 'process_low'
    
    conda "conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'quay.io/biocontainers/pigz:2.3.4' }"
	
    input:
        tuple val(meta), path(signatures) 
        
    output:
    
    path '*.csv'  emit: sig_csv
        
    
	script:
	"""
	
	"""
}