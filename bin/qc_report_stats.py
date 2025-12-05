#!/usr/bin/env python3

from genericpath import sameopenfile
from importlib.resources import path
import argparse
import pandas as pd

# Argument parser: get arguments from FAQCS and samtools text files
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample")
parser.add_argument("--stats", type=argparse.FileType("r"))
parser.add_argument("--base_content_before_trim", type=argparse.FileType("r"))
parser.add_argument("--base_content_after_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_before_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_after_trim", type=argparse.FileType("r"))
parser.add_argument("--reference", type=argparse.FileType("r"))
parser.add_argument("--samtools_coverage", type=argparse.FileType("r"))
parser.add_argument("--samtools_depth", type=argparse.FileType("r"))
parser.add_argument("--samtools_stats", type=argparse.FileType("r"))
parser.add_argument("--min_depth", type=int)
args = parser.parse_args()

# Sample name variable
sample_name = args.sample

# Reference length variable calculated from reference file
header = None
length = 0
for line in args.reference:
    # Trim newline
    line = line.rstrip()
    if line.startswith(">"):
        if header is not None:
            continue
        header = line[1:]
    else:
        length += len(line)

# Parse through stats.txt file for qc report variables
list = []
for lines in args.stats:
    list.append(lines)
reads_before_trim = list[1].split(" ")[2].strip("\n")
read_length_before_trim = list[3].split(" ")[2]
reads_after_trim = list[6].split(" ")[2].strip("\n")
reads_after_trim_percent = list[6].split(": ")[1].strip("\n")
read_length_after_trim = list[8].split(" ")[3]
paired_reads_after_trim = list[9].split(": ")[1].strip("\n")
unpaired_reads_after_trim = list[11].split(": ")[1].strip("\n")
coverage_numer_before = float(read_length_before_trim) * float(reads_before_trim)
coverage_before = coverage_numer_before / length
coverage_numer_after = float(read_length_after_trim) * float(reads_after_trim)
coverage_after = coverage_numer_after / length

# Calculate the GC content from base_conent.txt file using Pandas
df1 = pd.read_csv(
    args.base_content_before_trim, sep=r'\s+', header=None, index_col=None
)
df1_subset = df1[df1[0] == "GC"].copy()
df1_subset.loc[:, 3] = df1_subset[1] * df1_subset[2]
sum_reads = sum(df1_subset[2])
sum_reads_GC_content = sum(df1_subset[3])
# Formatting, 2 decimal places and adding percent sign
GC_content_before = (sum_reads_GC_content) / (sum_reads)
GC_content_before = "{:.2f}".format(GC_content_before)
GC_content_before = str(GC_content_before + "%")

# Calculate the GC content from base_conent.txt file using Pandas
df2 = pd.read_csv(
    args.base_content_after_trim, sep=r'\s+', header=None, index_col=None
)
df2_subset = df2[df2[0] == "GC"].copy()
df2_subset.loc[:, 3] = df2_subset[1] * df2_subset[2]
sum_reads = sum(df2_subset[2])
sum_reads_GC_content = sum(df2_subset[3])
# Formatting, 2 decimal places and adding percent sign
GC_content_after = (sum_reads_GC_content) / (sum_reads)
GC_content_after = "{:.2f}".format(GC_content_after)
GC_content_after = str(GC_content_after + "%")

# Calculate the average phred/quality score from qual_scores.txt file using Pandas
df3 = pd.read_csv(args.qual_scores_before_trim, sep=r'\s+', index_col=None)
df3["x"] = df3["Score"] * df3["readsNum"]
sum_reads_num = sum(df3["readsNum"])
phred_avg_before = sum(df3["x"]) / sum_reads_num
# Formatting. 2 decimal points
phred_avg_before = "{:.2f}".format(phred_avg_before)

# Calculate the average phred/quality score from qual_scores.txt file using Pandas
df4 = pd.read_csv(args.qual_scores_after_trim, sep=r'\s+', index_col=None)
df4["x"] = df4["Score"] * df4["readsNum"]
sum_reads_num = sum(df4["readsNum"])
phred_avg_after = sum(df4["x"]) / sum_reads_num
# Formatting. 2 decimal points
phred_avg_after = "{:.2f}".format(phred_avg_after)

# Parse samtools coverage output for mean depth and mapped reads
# Format: #rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
# The header line starts with # so we need to handle it specially
df_coverage = pd.read_csv(
    args.samtools_coverage,
    sep='\t',
    skiprows=1,  # Skip the header line that starts with #
    names=['rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq']
)

# Calculate weighted mean depth across all references
# Weight by the length of each reference (endpos - startpos)
df_coverage['ref_length'] = df_coverage['endpos'] - df_coverage['startpos']
total_ref_length = df_coverage['ref_length'].sum()
df_coverage['weighted_depth'] = df_coverage['meandepth'] * df_coverage['ref_length']
mean_depth_coverage = df_coverage['weighted_depth'].sum() / total_ref_length

# Parse samtools stats for total reads and reads mapped
# Format: SN	raw total sequences:	<number>
# Format: SN	reads mapped:	<number>
total_reads = 0
reads_mapped_count = 0
for line in args.samtools_stats:
    if line.startswith('SN') and 'raw total sequences:' in line:
        # Line format: SN\traw total sequences:\t<number>
        parts = line.strip().split('\t')
        total_reads = int(parts[2])
    elif line.startswith('SN') and 'reads mapped:' in line:
        # Line format: SN\treads mapped:\t<number>
        parts = line.strip().split('\t')
        reads_mapped_count = int(parts[2])

# Calculate percentage of reads mapped using total reads from BAM file
if total_reads > 0:
    reads_mapped_percentage = (reads_mapped_count / total_reads) * 100
else:
    reads_mapped_percentage = 0.0
reads_mapped = f"{reads_mapped_count} ({reads_mapped_percentage:.2f}%)"

# Parse samtools depth output to calculate genome fraction at min_depth
# Format: reference	position	depth
df_depth = pd.read_csv(args.samtools_depth, sep='\t', header=None, names=['reference', 'position', 'depth'])

# Count positions with depth >= min_depth
positions_above_threshold = len(df_depth[df_depth['depth'] >= args.min_depth])
total_positions = length  # Use reference length calculated earlier

# Calculate percentage
if total_positions > 0:
    genome_fraction_value = (positions_above_threshold / total_positions) * 100
    genome_fraction = "{:.2f}".format(genome_fraction_value) + "%"
else:
    genome_fraction = "0.00%"

# Preparing output list with variables and then reformatting into a string
output_list = [
    sample_name,
    str(reads_before_trim),
    str(GC_content_before),
    str(phred_avg_before),
    "{:.2f}".format(coverage_before),
    str(reads_after_trim_percent),
    str(paired_reads_after_trim),
    str(unpaired_reads_after_trim),
    str(GC_content_after),
    str(phred_avg_after),
    "{:.2f}".format(coverage_after),
    "{:.2f}".format(mean_depth_coverage),
    str(reads_mapped),
    genome_fraction
]

# Creating tab delimited string for qc report generation
print('\t'.join(output_list))
