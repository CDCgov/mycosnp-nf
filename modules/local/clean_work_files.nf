/**
 * Clean intermediate files.
 Taken from: https://github.com/SystemsGenetics/GEMmaker
 */
process clean_work_files {
    tag { meta.id }
    label "local"

    input:
    tuple val(meta), val(files)

    output:
    val(1), emit: IS_CLEAN

    script:
    """
    clean_work_files.sh "${files.join(" ")}"
    """
}