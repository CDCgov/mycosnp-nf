process FASTA_MERGE {
    
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    path("combined/alignment.fasta"), emit: alignment

    script:
    def use_zcat = false
    def zcat_files = ""
    def use_cat = false
    def cat_files = ""

    for (int i=0; i < fasta.size(); i++)
    {
        thisFasta = fasta[i]
        if (thisFasta.getName().endsWith(".gz")) 
        {
            use_zcat = true
            zcat_files += " $thisFasta "
        } else
        {
            use_cat = true
            cat_files += " $thisFasta "
        }
    }

    """
    echo ${meta.id} $cat_files
    echo ${meta.id} $zcat_files
    mkdir combined

    if [[ "$use_cat" == "true" ]]; then
        cat $cat_files >> combined/alignment.fasta
    fi
    if [[ "$use_zcat" == "true" ]]; then
        zcat $zcat_files >> combined/alignment.fasta
    fi
    """

}
