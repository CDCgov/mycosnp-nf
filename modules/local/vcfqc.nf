process VCF_QC {
    tag "vcf-qc"
    label 'process_low'
 
    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(vcffasta)

    output:
    path("vcf-qc-report.txt"),   emit: vcf_qc_report
	

	script:
	"""
    printf "Sample Name\\tLength\\tNumber-N\\n" > vcf-qc-report.txt
    awk '\$0 ~ ">" {if (NR > 1) {print c "\\t" d;} c=0;d=0;printf substr(\$0,2,200) "\\t"; } \$0 !~ ">" {c+=length(\$0);d+=gsub(/N/, "");d+=gsub(/n/, "")} END { print c "\\t" d; }' $vcffasta >> vcf-qc-report.txt
	"""
}
