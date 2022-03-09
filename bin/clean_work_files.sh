#!/usr/bin/env bash
# This script is meant for cleaning any file in a Nextflow work directory.
# The $files_list variable is set within the Nextflow process and should
# contain the list of files that need cleaning. This can be done by creating
# a channel in a process that creates files, and merging that channel with
# a signal from another process indicating the files are ready for cleaning.
#
# The cleaning process empties the file, converts it to a sparse file so it
# has an acutal size of zero but appears as the original size, the access
# and modify times are kept the same.
#
# Taken from: https://github.com/SystemsGenetics/GEMmaker
files_list="$1"

for file in ${files_list}; do
  # Remove cruff added by Nextflow
  file=`echo $file | perl -p -e 's/[\\[,\\]]//g'`
  if [ -e $file ]; then
    # Log some info about the file for debugging purposes
    echo "cleaning $file"
    stat $file
    # Get file info: size, access and modify times
    size=`stat --printf="%s" $file`
    atime=`stat --printf="%X" $file`
    mtime=`stat --printf="%Y" $file`

    # Make the file size 0 and set as a sparse file
    > $file
    truncate -s $size $file
    # Reset the timestamps on the file
    touch -a -d @$atime $file
    touch -m -d @$mtime $file
  fi
done