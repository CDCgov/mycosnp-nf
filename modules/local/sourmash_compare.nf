process SOURMASH_COMPARE {
    tag "$meta.id"
    label 'process_low'
    //copied from nf-core modules and modified

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.4--hdfd78af_0':
        'biocontainers/sourmash:4.8.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(signatures)
    path file_list // optional file
    val save_numpy_matrix
    val save_csv

    output:
    tuple path("*comp.npy.labels.txt"), path("*comp.npy"), optional:true, emit: matrix
    tuple val(meta), path("*jaccard_comp.csv")       , optional:true, emit: jaccard_csv //for all v all
    tuple val(meta), path("*ani_comp.csv")           , optional:true, emit: ani_csv //for all v all
    
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args      ?: '--ksize 31'
    def ani_args       = '--ksize 31 --ani'
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def comp           = save_numpy_matrix  ? "--output comp.npy"  : ''
    def jaccard_csv    = save_csv           ? "--csv jaccard_comp.csv" : ''
    def ani_csv        = save_csv           ? "--csv ani_comp.csv" : ''
    if ( !save_numpy_matrix && !save_csv ) error "Supply either save_numpy_matrix, save_csv, or both or no output will be created"
    def ffile = file_list //? "--from-file ${file_list}" : ''
    def sigs = signatures ? "${signatures.join(' ')}" : ''
    if ( !file_list && !signatures ) error "Supply either signatures, file_list, or both"
    """
    sourmash compare \\
        $ani_args \\
        --processes ${task.cpus} \\
        ${comp} \\
        --csv ${prefix}.ani_comp.csv \\
        ${sigs} \\
        ${ffile} \\
        
    
    sourmash compare \\
        $args \\
        --processes ${task.cpus} \\
        ${comp} \\
        --csv ${prefix}.jacc_comp.csv \\
        ${sigs} \\
        ${ffile} \\
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}