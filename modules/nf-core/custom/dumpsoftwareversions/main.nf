process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda (params.enable_conda ? 'bioconda::multiqc=1.27.1' : null)
    container 'quay.io/staphb/multiqc:1.27.1'

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
