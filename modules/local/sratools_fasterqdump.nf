process SRATOOLS_FASTERQDUMP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/sratools.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.0.8--h9f5acd7_0' :
        'biocontainers/sra-tools:3.0.8--h9f5acd7_0' }"

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path(fastq_output), emit: reads
    tuple val(meta), path(md5_output)  , emit: md5
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    // Paired-end data extracted by fasterq-dump (--split-3 the default) always creates
    // *_1.fastq *_2.fastq files but sometimes also an additional *.fastq file
    // for unpaired reads which we ignore here.
    fastq_output = meta.single_end ? '*.fastq.gz'     : '*_{1,2}.fastq.gz'
    md5_output   = meta.single_end ? '*.fastq.gz.md5' : '*_{1,2}.fastq.gz.md5'
    """
    eval "\$(vdb-config -o n NCBI_SETTINGS | sed 's/[" ]//g')"
    if [[ ! -f "\${NCBI_SETTINGS}" ]]; then
        mkdir -p "\$(dirname "\${NCBI_SETTINGS}")"
        printf '${config}' > "\${NCBI_SETTINGS}"
    fi
    fasterq-dump \\
        $args \\
        --threads $task.cpus \\
        ${sra.name}
    pigz \\
        $args2 \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq
    ## Rename FastQ files by meta.id
    if [ -f  ${sra.name}.fastq.gz ]; then
        mv ${sra.name}.fastq.gz ${meta.id}.fastq.gz
        md5sum ${meta.id}.fastq.gz > ${meta.id}.fastq.gz.md5
    fi
    if [ -f  ${sra.name}_1.fastq.gz ]; then
        mv ${sra.name}_1.fastq.gz ${meta.id}_1.fastq.gz
        md5sum ${meta.id}_1.fastq.gz > ${meta.id}_1.fastq.gz.md5
    fi
    if [ -f  ${sra.name}_2.fastq.gz ]; then
        mv ${sra.name}_2.fastq.gz ${meta.id}_2.fastq.gz
        md5sum ${meta.id}_2.fastq.gz > ${meta.id}_2.fastq.gz.md5
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
// taken from [nf-core/fetchngs](https://github.com/nf-core/fetchngs/blob/master/modules/local/sratools_fasterqdump.nf)
