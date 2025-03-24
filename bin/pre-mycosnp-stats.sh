#!/bin/bash

# pre-mycosnp-stats.sh v1.1.0
# Authors:
# Jared Johnson, jared.johnson@doh.wa.gov
# Zack Mudge, ZMudge@cdc.gov

version="v1.1.0"

#----- HELP & VERSION -----#
# help message
if [[ "$1" == "-h" || "$1" == "--help" || "$1" == "-help" ]]; then
    echo -e "pre-mycosnp-stats.sh [-r path/to/ref_assembly (optional)] [sample_name] [path/to/sample_assembly] [path/to/faqcs/*.stats.txt] [path/to/faqcs/*.for_qual_histogram.txt]" && exit 0
fi

# version
if [[ "$1" == "-v" || "$1" == "--version" || "$1" == "-version" ]]; then
    echo -e ${version} && exit 0
fi

# Initialize variables
sample=""
assembly=""
ref=""
faqcs_stats=""
faqcs_qual=""

# Parse command line options using getopts
while getopts ":r:" opt; do
    case $opt in
        r)
            ref="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done

shift $((OPTIND - 1)) # Shift off the processed options

# Assign remaining positional arguments after options are processed.
if [[ $# -lt 4 ]]; then
    echo "Error: Missing required arguments."
    exit 1
fi

sample="$1"
assembly="$2"
faqcs_stats="$3"
faqcs_qual="$4"

# decompress the sample and/or reference assemblies
gzip -d "${assembly}" "${ref}" || true

#----- READ STATS -----#
# total trimmed reads
trmd_reads=$(cat "${faqcs_stats}" | grep 'Reads #:' | sed -n 2p | cut -f 3 -d ' ')
# average trimmed read Phred score, to 2 decimal places
avg_phred=$(cat "${faqcs_qual}" | awk '{print $2,$1*$2}' | awk '{reads += $1} {qual += $2} END {printf "%.2f", qual/reads}')

#----- ASSEMBLY STATS -----#
# sample assembly length
sample_length=$(cat "${assembly%.gz}" | grep -v ">" | tr -d '\n\t\r ' | wc -c)

#----- REF STATS -----#
if [[ -n "$ref" ]]; then # Check if ref is provided (not empty)
    # reference assembly length
    ref_length=$(cat "${ref%.gz}" | grep -v ">" | tr -d '\n\t\r ' | wc -c)

    # reference %GC, to 2 decimal places
    ref_gc_count=$(cat "${ref%.gz}" | grep -v ">" | grep -Eo "G|C" | wc -l)
    ref_gc_perc=$(echo -e "${ref_gc_count}\t${ref_length}" | awk '{printf "%.2f", 100*$1/$2}')

    # estimated average depth of coverage, to 1 decimal place
    est_depth=$(echo "$trmd_reads $ref_length" | awk '{printf "%.1f", $1 / $2}')
else
    # Set variables to empty strings if ref is not provided.
    ref_length=""
    est_depth=""
    ref_gc_perc=""
fi

# sample assembly % GC, to 2 decimal places
sample_gc_count=$(cat "${assembly%.gz}" | grep -v ">" | grep -Eo "G|C" | wc -l)
sample_gc_perc=$(echo -e "${sample_gc_count}\t${sample_length}" | awk '{printf "%.2f", 100*$1/$2}')

#----- OUTPUT -----#
echo "${trmd_reads},${avg_phred},${sample_length},${sample_gc_perc},${ref_length},${est_depth},${ref_gc_perc}"
