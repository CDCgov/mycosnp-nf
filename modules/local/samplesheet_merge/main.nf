process SAMPLESHEET_MERGE {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::perl=5.22.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.22.2.1' :
        'quay.io/biocontainers/perl:5.22.2.1' }"

    input:
    path(samplesheet)

    output:
    path('*.system.csv'), emit: csv
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in bin/
    """
    mycosnp_combine_lanes.pl -i ${samplesheet} > samplesheet.system.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v 2>&1 | grep -o 'v[0-9.][0-9.]*' | sed 's/v//')
    END_VERSIONS
    """

    stub:
    """
    touch samplesheet.system.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v 2>&1 | grep -o 'v[0-9.][0-9.]*' | sed 's/v//')
    END_VERSIONS
    """
}
