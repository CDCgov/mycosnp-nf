process SNPEFFR {
	
	container "ghcr.io/cdcgov/snpeffr:master"
	
	input:
	tuple val(meta), path(input)

	output:
	path '*.csv'		,emit: report
	path "versions.yml"	,emit: versions

	when:
	task.ext.when == null || task.ext.when    

	script:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"

	"""
	Rscript $PWD/mycosnp-nf-1.5/bin/snpeffr.R -f $input \\
	$args -o ${prefix}.csv
	
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	snpeffr: \$(snpeffr --version | sed 's/v//')
	END_VERSIONS
	"""
}
