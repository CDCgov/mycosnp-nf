#!/bin/bash

species=$1
db=$2
assembly=$3
prefix=$4

# select a mash sketch file (if listed)
sketch_file=$(cat ${db}/sketch_list.csv | tr ',' '\t' | awk -v s=${species} '$1 == s {print $2}')

# run mash if sketch exists
if [[ -f "${db}/${sketch_file}" ]]
then
    echo -e "sample,subtype,mash_dist,est_ANI" > ${prefix}_subtype.csv
    mash dist ${db}/${sketch_file} ${assembly} | sort -nk 3 | awk -v s=${prefix} 'BEGIN {OFS = ","} {print s,$1,$3,100*(1-$3)}' >> ${prefix}_subtype.csv
else
    touch ${prefix}_subtype.csv
fi
