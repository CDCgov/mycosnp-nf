process GET_QC_REF {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ncbi-datasets=14.26.0" : null)
    container 'staphb/ncbi-datasets:15.2.0'

    input:
    tuple val(meta), path(gambit_results)

    output:
    tuple val(meta), path("*.fna"), emit: qc_ref

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # determine species call from Gambit
    rank=\$(cat ${gambit_results} | tr ',' '\t' | grep -v "predicted.rank" | cut -f 3)
    accession=\$(cat ${gambit_results} | tr ',' '\t' | grep -v 'closest.description' | cut -f 7 | cut -f 1 -d ' ' | tr -d '[] \t\n\r')

    if [[ "\${rank}" == "species" ]]
    then
        echo "datasets download genome accession \${accession}"
        datasets download genome accession \${accession} && unzip ncbi_dataset.zip && mv ncbi_dataset/data/*/*.fna ./
    fi
    """
}
