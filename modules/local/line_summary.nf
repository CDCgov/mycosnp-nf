process LINE_SUMMARY {
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::ncbi-datasets=14.26.0" : null)
    container 'staphb/ncbi-datasets:15.2.0'

    input:
    tuple val(meta), path(assembly), path(faqcs), path(gambit), path(subtype)

    output:
    tuple val(meta), path("*_linesummary.csv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "$meta.id"
    """
    ## Extract relevant fields from Gambit output
    taxa=\$(cat ${gambit} | tr ',' '\t' | grep -v "predicted.name" | cut -f 2)
    rank=\$(cat ${gambit} | tr ',' '\t' | grep -v "predicted.rank" | cut -f 3)
    accession=\$(cat ${gambit} | tr ',' '\t' | grep -v 'closest.description' | cut -f 7 | cut -f 1 -d ' ' | tr -d '[] \t\n\r')
    
    # Download the Gambit reference for estimating average coverage
    if [[ "\${rank}" == "species" ]]
    then
        echo "datasets download genome accession \${accession}"
        datasets download genome accession \${accession} && unzip ncbi_dataset.zip && mv ncbi_dataset/data/*/*.fna ./ref.fa
    else
        echo "Error: Species could not be assigned - check the Gambit output for more details" && exit 1
    fi

    # Gather QC stats
    fast-stats.sh \\
        ${prefix} \\
        ${assembly} \\
        ref.fa \\
        *.stats.txt \\
        *.for_qual_histogram.txt \\
        > stats_cols

    # Extract subtype info if available
    if [[ -f ${subtype} ]]
    then
        subtype_call=\$(head -n 2 ${subtype} | grep -v 'sample,subtype,mash_dist,est_ANI' | tr ',' '\t' | cut -f 2)
        subtype_ani=\$(head -n 2 ${subtype} | grep -v 'sample,subtype,mash_dist,est_ANI' | tr ',' '\t' | cut -f 4)
    else
        subtype_call="NA"
        subtype_ani="NA"
    fi

    # Create line summary
    echo "\${taxa},\${subtype_call},\${subtype_ani},\${accession}" > taxa_cols
    paste -d ',' stats_cols taxa_cols > ${prefix}_linesummary.csv
    
    """
}
