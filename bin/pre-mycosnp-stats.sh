#!/bin/bash

# pre-mycosnp-stats.sh v1.0.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

version="v1.0.0"

# inputs
sample=$1
assembly=$2
ref=$3
faqcs_stats=$4
faqcs_qual=$5

#----- HELP & VERSION -----#
# help message
if [ ${sample} == "-h" ] || [ ${sample} == "--help" ] || [ ${sample} == "-help" ]
then
    echo -e "subtyper.sh [sample_name] [path/to/sample_assembly] [path/to/ref_assembly] [path/to/faqcs/*.stats.txt] [path/to/faqcs/*.for_qual_histogram.txt]" && exit 0
fi

# version
if [ ${sample} == "-v" ] || [ ${sample} == "--version" ] || [ ${sample} == "-version" ]
then
    echo -e ${version} && exit 0
fi

# decompress the sample and/or reference assemblies - if needed
gzip -d ${assembly} ${ref} || true

#----- ASSEMBLY STATS -----#
# sample assembly length
sample_length=$(cat ${assembly%.gz} | grep -v ">" | tr -d '\n\t\r ' | wc -c)
# reference assembly length
ref_length=$(cat ${ref%.gz} | grep -v ">" | tr -d '\n\t\r ' | wc -c)
# sample % GC
sample_gc_count=$(cat ${assembly%.gz} | grep -v ">" | grep -Eo "G|C" | wc -l)
sample_gc_perc=$(echo -e "${sample_gc_count}\t${sample_length}" | awk '{print 100*$1/$2}')
# reference %GC
ref_gc_count=$(cat ${ref%.gz} | grep -v ">" | grep -Eo "G|C" | wc -l)
ref_gc_perc=$(echo -e "${ref_gc_count}\t${ref_length}" | awk '{print 100*$1/$2}')

#----- READ STATS -----#
# total trimmed reads
trmd_reads=$(cat ${faqcs_stats} | grep 'Reads #:' | sed -n 2p | cut -f 3 -d ' ')
# total trimmed bases
trmd_bases=$(cat ${faqcs_stats} | grep 'Total bases:' | sed -n 2p | cut -f 3 -d ' ')
# average trimmed read Phred score
avg_phred=$(cat ${faqcs_qual} | awk '{print $3,$1*$3}' | awk '{bases += $1} {qual += $2} END {print qual/bases}')
# estimated average depth of coverage
est_depth=$((trmd_bases/ref_length))

#----- OUTPUT -----#
echo "${sample},${trmd_reads},${avg_phred},${est_depth},${sample_length},${ref_length},${sample_gc_perc},${ref_gc_perc}"
