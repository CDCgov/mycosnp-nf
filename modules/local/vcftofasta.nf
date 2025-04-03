process VCF_TO_FASTA {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::scipy=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scipy%3A1.1.0' :
        'quay.io/biocontainers/scipy:1.1.0' }"

    input:
    tuple val(meta), path(vcf), path(samplelist), val(max_amb_samples), val(max_perc_amb_samples), val(min_depth)
    path(fasta)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    // path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def is_compressed_vcf = vcf.getName().endsWith(".gz") ? true : false
    def vcf_name = vcf.getName().replace(".gz", "")

    """
    NUM_SAMPLES=\$(cat $samplelist | wc -l)
    if [ $max_perc_amb_samples > 0 ]; then
        MAX_AMB_SAMPLES=\$(echo "\${NUM_SAMPLES} $max_perc_amb_samples" | awk '{x=\$1*(\$2/100); y=int(x); x=(y<1?1:y)} END {print x}')
    else
        MAX_AMB_SAMPLES=$max_amb_samples
    fi

    if [ "$is_compressed_vcf" == "true" ]; then
        gzip -c -d $vcf > $vcf_name
    fi

    vcfSnpsToFasta.py --max_amb_samples \$MAX_AMB_SAMPLES --min_depth $min_depth $vcf_name > ${prefix}_vcf-to-fasta.fasta
    echo "NUM_SAMPLES=\$NUM_SAMPLES" >> log.txt
    echo "MAX_PERC_AMB_SAMPLES=$max_perc_amb_samples" >> log.txt
    echo "MAX_AMB_SAMPLES=\$MAX_AMB_SAMPLES" >> log.txt
    echo "MIN_DEPTH=$min_depth" >> log.txt
    
    """

}
