#!/bin/bash

# inputs
sample=$1
assembly=$2
ref=$3
faqcs_stats=$4
faqcs_qual=$5

# decompress any input files
gzip -d *.gz || true

# read stats
trmd_reads=$(cat ${faqcs_stats} | grep 'Reads #:' | sed -n 2p | cut -f 3 -d ' ')
trmd_bases=$(cat ${faqcs_stats} | grep 'Total bases:' | sed -n 2p | cut -f 3 -d ' ')

avg_phred=$(cat ${faqcs_qual} | awk '{print $3,$1*$3}' | awk '{bases += $1} {qual += $2} END {print qual/bases}')

# assembly stats
sample_length=$(cat ${assembly%.gz} | grep -v ">" | tr -d '\n\t\r ' | wc -c)
ref_length=$(cat ${ref%.gz} | grep -v ">" | tr -d '\n\t\r ' | wc -c)

sample_gc_count=$(cat ${assembly%.gz} | grep -v ">" | grep -Eo "G|C" | wc -l)
ref_gc_count=$(cat ${ref%.gz} | grep -v ">" | grep -Eo "G|C" | wc -l)

sample_gc_perc=$(echo -e "${sample_gc_count}\t${sample_length}" | awk '{print 100*$1/$2}')
ref_gc_perc=$(echo -e "${ref_gc_count}\t${ref_length}" | awk '{print 100*$1/$2}')

est_depth=$((trmd_bases/ref_length))

echo "${sample},${trmd_reads},${avg_phred},${est_depth},${sample_length},${ref_length},${sample_gc_perc},${ref_gc_perc}"
