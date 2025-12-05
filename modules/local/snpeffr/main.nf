process SNPEFFR {
	tag "$meta.id"
	label 'process_medium'

	container "ghcr.io/cdcgov/snpeffr:v1.1.1"

	input:
	tuple val(meta), path(input)
	val(ref_name)

	output:
	path('*.csv'), 		 emit: report
	path "versions.yml", emit: versions

	when:
	task.ext.when == null || task.ext.when

	script:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"
	"""
	snpeffr.R -f ${input} \\
	${args} -o ${prefix}_${ref_name}.csv

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	snpeffr: \$(snpeffr --version | sed 's/v//')
	END_VERSIONS
	"""

	stub:
	def prefix = task.ext.prefix ?: "${meta.id}"
	"""
	touch ${prefix}_${ref_name}.csv

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	snpeffr: \$(snpeffr --version | sed 's/v//')
	END_VERSIONS
	"""
}
