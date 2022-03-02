process VCF_QC {
    tag "vcf-qc"
    label 'process_low'
 
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
