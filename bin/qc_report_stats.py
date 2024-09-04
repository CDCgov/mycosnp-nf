#!/usr/bin/env python3
from genericpath import sameopenfile
from importlib.resources import path
import argparse
import pandas as pd

# Arguement parser: get arguments from FAQCS text files
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample")
parser.add_argument("--stats", type=argparse.FileType("r"))
parser.add_argument("--base_content_before_trim", type=argparse.FileType("r"))
parser.add_argument("--base_content_after_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_before_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_after_trim", type=argparse.FileType("r"))
parser.add_argument("--reference", type=argparse.FileType("r"))
parser.add_argument("--bam_coverage", type=argparse.FileType("r"))
parser.add_argument("--genome_fraction", type=argparse.FileType("r"))
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
    args.base_content_before_trim, delim_whitespace=True, header=None, index_col=None
)
df1_subset = df1[df1[0] == "GC"]
df1_subset[3] = df1[1] * df1[2]
sum_reads = sum(df1_subset[2])
sum_reads_GC_content = sum(df1_subset[3])
# Formatting, 2 decimal places and adding percent sign
GC_content_before = (sum_reads_GC_content) / (sum_reads)
GC_content_before = "{:.2f}".format(GC_content_before)
GC_content_before = str(GC_content_before + "%")

# Calculate the GC content from base_conent.txt file using Pandas
df2 = pd.read_csv(
    args.base_content_after_trim, delim_whitespace=True, header=None, index_col=None
)
df2_subset = df2[df2[0] == "GC"]
df2_subset[3] = df2[1] * df2[2]
sum_reads = sum(df2_subset[2])
sum_reads_GC_content = sum(df2_subset[3])
# Formatting, 2 decimal places and adding percent sign
GC_content_after = (sum_reads_GC_content) / (sum_reads)
GC_content_after = "{:.2f}".format(GC_content_after)
GC_content_after = str(GC_content_after + "%")

# Calculate the average phred/quality score from qual_scores.txt file using Pandas
df3 = pd.read_csv(args.qual_scores_before_trim, delim_whitespace=True, index_col=None)
df3["x"] = df3["Score"] * df3["readsNum"]
sum_reads_num = sum(df3["readsNum"])
phred_avg_before = sum(df3["x"]) / sum_reads_num
# Formatting. 2 decimal points
phred_avg_before = "{:.2f}".format(phred_avg_before)

# Calculate the average phred/quality score from qual_scores.txt file using Pandas
df4 = pd.read_csv(args.qual_scores_after_trim, delim_whitespace=True, index_col=None)
df4["x"] = df4["Score"] * df4["readsNum"]
sum_reads_num = sum(df4["readsNum"])
phred_avg_after = sum(df4["x"]) / sum_reads_num
# Formatting. 2 decimal points
phred_avg_after = "{:.2f}".format(phred_avg_after)

# Parsing coverage (mean depth) and percent coverage of reference sequence from Qualimap bamqc report genome_results.txt file
for line in args.bam_coverage:
    if line.__contains__("mean coverageData"):
        mean_depth_coverage = float(line.split("= ")[1].strip("X\n"))
    if line.__contains__("number of mapped reads"):
        reads_mapped = line.split("= ")[1].replace(",","").strip("\n")

# Parsing genome fraction from qualimap genome_fraction_coverage.txt
next(args.genome_fraction)
for line in args.genome_fraction:
    columns = line.strip().split('\t')
    if float(columns[0]) == float(args.min_depth):
        genome_fraction = float(columns[1])
        break
    else:
        genome_fraction = 0

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
    "{:.2f}".format(genome_fraction)
]

# Creating tab delimited string for qc report generation
print('\t'.join(output_list))
