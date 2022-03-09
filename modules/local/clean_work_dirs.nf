/**
 * Clean intermediate directories.
  Taken from: https://github.com/SystemsGenetics/GEMmaker
 */
process clean_work_dirs {
    tag { meta.id }
    label "local"

    input:
    tuple val(meta), val(directory)

    output:
    val(1), emit: IS_CLEAN

    script:
    """
    for dir in ${directory}; do
      if [ -e \$dir ]; then
        echo "Cleaning: \$dir"
        files=`find \$dir -type  f `
        echo "Files to delete: \$files"
        clean_work_files.sh "\$files" "null"
      fi
    done
    """
}