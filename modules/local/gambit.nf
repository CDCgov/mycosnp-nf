process GAMBIT_QUERY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gambit=1.1.0=py312h0fa9677_2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gambit:1.1.0--py312h0fa9677_2' :
        'quay.io/biocontainers/gambit:1.1.0--py312h0fa9677_2' }"
    
    input:
    tuple val(meta), path(assembly)
    path db_file
    path h5_files_dir

    output:
    tuple val(meta), path("*_gambit.csv"), emit: taxa
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_h5_file = "gambit_signatures.gs"

    """
    # Re-combine chunks of GAMBIT signature file
    > ${output_h5_file}  # Create or clear the combined signatures file
    for file in ${h5_files_dir}/*; do
        if [[ "\$file" == *.gz ]]; then
            gunzip -c "\$file" >> ${output_h5_file}
        else
            cat "\$file" >> ${output_h5_file}
        fi
    done

    gambit ${args} -d ./ query ${assembly} > ${prefix}_gambit.csv

    rm ${output_h5_file}
    # Replace "Candidozyma auris" with "Candida auris" in the output CSV
    sed -i 's/Candidozyma auris/Candida auris/g' ${prefix}_gambit.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gambit: \$(gambit --version | tr -d "[a-z, ]")
    END_VERSIONS
    """
}
