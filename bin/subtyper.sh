#!/bin/bash

# subtyper.sh v1.0.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

version="v1.0.0"

# inputs
sample=$1
species=${2// /_}
db=${3%/}
seq=$4

#----- HELP & VERSION -----#
# help message
if [ ${sample} == "-h" ] || [ ${sample} == "--help" ] || [ ${sample} == "-help" ]
then
    echo -e "subtyper.sh [sample_name] [species] [path/to/db] [path/to/assembly|reads]" && exit 0
fi

# version
if [ ${sample} == "-v" ] || [ ${sample} == "--version" ] || [ ${sample} == "-version" ]
then
    echo -e ${version} && exit 0
fi

#----- SUBTYPER -----#
# select a mash sketch file (if listed)
if [ -f "${db}/sketch_list.csv" ]
then
    sketch_file=$(cat ${db}/sketch_list.csv | tr ',' '\t' | awk -v s=${species} '$1 == s {print $2}')
else
    echo "Error: ${db}/sketch_list.csv does not exist." && exit 1
fi

# run mash if sketch exists
if [ -f "${db}/${sketch_file}" ]
then
    echo -e "sample,subtype,mash_dist,est_ANI" > ${name}_subtype.csv
    mash dist ${db}/${sketch_file} ${seq} | sort -nk 3 | awk -v s=${sample} 'BEGIN {OFS = ","} NR==1 {print s,$1,$3,100*(1-$3)}' >> ${sample}_subtype.csv
else
    echo -e "Error: Sketch file for ${species} not found."
fi
