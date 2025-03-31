process PRE_MYCOSNP_INDV_SUMMARY {
    tag "$meta.id"
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
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ## Extract relevant fields from Gambit output
    taxon=\$(cat "${gambit}" | grep -v "predicted.name" | cut -f 2 -d ',')
    rank=\$(cat "${gambit}" | grep -v "predicted.rank" | cut -f 3 -d ',')
    distance=\$(printf "%.4f" \$(cat "${gambit}" | grep -v "closest.distance" | cut -f 6 -d ','))
    closest=\$(cat "${gambit}" | grep -v 'closest.description' | cut -f 7 -d ',')
    closest_accession=\$(echo "\${closest}" | cut -f 1 -d ' ' | tr -d '[] \\t\\n\\r')

    subtype_call=""
    subtype_accession=""
    subtype_ani=""

    if [[ "\${rank}" == "species" ]]; then
        # Download the Gambit reference for estimating average depth of coverage
        retry_with_backoff.sh -d 15 datasets download genome accession \${closest_accession}
        unzip ncbi_dataset.zip && mv ncbi_dataset/data/*/*.fna ./ref.fa

        # Gather QC stats
        pre-mycosnp-stats.sh \\
            -r ref.fa \\
            "${prefix}" \\
            "${assembly}" \\
            "${prefix}.stats.txt" \\
            "${prefix}.for_qual_histogram.txt" \\
            > stats_cols

        # Extract subtype info if available
        if [ -s "${prefix}_subtype.csv" ]; then
            second_col=\$(head -n 2 "${subtype}" | grep -v 'sample,subtype,mash_dist,est_ANI' | cut -f 2 -d ',')
            subtype_call=\$(echo "\${second_col}" | awk -F'-' '{print substr(\$1, 2)}')
            subtype_accession=\$(echo "\${second_col}" | awk -F'-' '{print \$2}')
            subtype_ani=\$(head -n 2 "${subtype}" | grep -v 'sample,subtype,mash_dist,est_ANI' | cut -f 4 -d ',')

            # Check if taxon is "Candida auris" and if subtype_ani is less than 99.7 using awk
            if [[ "\${taxon}" == "Candida auris" ]] && awk -v ani="\${subtype_ani}" 'BEGIN { exit !(ani < 99.7) }'; then
                subtype_call="ANI is less than the established Candida auris clade separation threshold of 99.7"
            fi
        fi
    elif [[ "\${rank}" == "genus" || "\${rank}" == "" ]]; then
        # Gather QC stats
        pre-mycosnp-stats.sh \\
            "${prefix}" \\
            "${assembly}" \\
            "${prefix}.stats.txt" \\
            "${prefix}.for_qual_histogram.txt" \\
            > stats_cols
    fi

    reported_rank="\${rank}"
    if [[ "\${rank}" == "" ]]; then
        reported_rank="no prediction"
    fi

    # Create line summary. Wrap `closest` in double quotes in case value contains commas.
    echo "${prefix},\${reported_rank},\${taxon},\${subtype_call},\${subtype_accession},\${subtype_ani},\${closest},\${distance}" > taxon_cols
    paste -d ',' taxon_cols stats_cols > "${prefix}_linesummary.csv"
    """
}
