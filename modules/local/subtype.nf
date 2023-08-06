process SUBTYPE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mash=2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash%3A2.3--he348c14_1' : 
        'quay.io/staphb/mash:2.3' }"
        
    input:
    tuple val(meta), path(gambit_results), path(assembly)
    path subtype_db

    output:
    path("*_subtype.txt"), emit: subtype

    shell:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    '''
    # determine species call from Gambit
    species=$(cat !{gambit_results} | tr ',' '\t' | cut -f 2)
    
    # select a mash sketch file (if listed)
    sketch_file=$(cat sketch_list.csv | tr ',' '\t' | awk -v s=${species} '$1 == s {print $2}')

    if [[ ${sketch_file} != "" ]]
    then
        echo -e "clade\tmash_dist\test_ANI" > !{prefix}_subtype.txt
        mash dist ${sketch_file} !{assembly} | sort -nk 3 | awk '{print $1"\t"$3"\t"100*(1-$3)}' | sed 's/.f*$//g' >> !{prefix}_subtype.txt
    else
        touch !{prefix}_subtype.txt
    fi
    '''
}
