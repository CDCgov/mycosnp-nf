#!/usr/bin/env python
from genericpath import sameopenfile
from importlib.resources import path
import argparse
import pandas as pd

# Arguement parser
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample")
parser.add_argument("--stats", type=argparse.FileType("r"))
parser.add_argument("--base_content_before_trim", type=argparse.FileType("r"))
parser.add_argument("--base_content_after_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_before_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_after_trim", type=argparse.FileType("r"))
parser.add_argument("--reference", type=argparse.FileType("r"))
args = parser.parse_args()

# Sample name variable
sample_name = args.sample

# Reference length variable
header = None
length = 0
for line in args.reference:
    # Trim newline
    line = line.rstrip()
    if line.startswith(">"):
        # If we captured one before, print it now
        if header is not None:
            print(header, length)
            length = 0
        header = line[1:]
    else:
        length += len(line)

# Stats report variables
list = []
for lines in args.stats:
    list.append(lines)
reads_before_trim = list[1].split(" ")[2].strip("\n")
read_length_before_trim = list[3].split(" ")[2]
reads_after_trim = list[6].split(": ")[1].strip("\n")
read_length_after_trim = list[8].split(" ")[3]
paired_reads_after_trim = list[9].split(": ")[1].strip("\n")
unpaired_reads_after_trim = list[11].split(": ")[1].strip("\n")
coverage_numer = float(read_length_before_trim) * float(reads_before_trim)
coverage = coverage_numer / length

# base_content_before_trim
df1 = pd.read_csv(
    args.base_content_before_trim, delim_whitespace=True, header=None, index_col=None
)
df1_subset = df1[df1[0] == "GC"]
df1_subset[3] = df1[1] * df1[2]
sum_reads = sum(df1_subset[2])
sum_reads_GC_content = sum(df1_subset[3])
GC_content_before = (sum_reads_GC_content) / (sum_reads)
GC_content_before = "{:.2f}".format(GC_content_before)
GC_content_before = str(GC_content_before + "%")

# base_content_after_trim
df2 = pd.read_csv(
    args.base_content_after_trim, delim_whitespace=True, header=None, index_col=None
)
df2_subset = df2[df2[0] == "GC"]
df2_subset[3] = df2[1] * df2[2]
sum_reads = sum(df2_subset[2])
sum_reads_GC_content = sum(df2_subset[3])
GC_content_after = (sum_reads_GC_content) / (sum_reads)
GC_content_after = "{:.2f}".format(GC_content_after)
GC_content_after = str(GC_content_after + "%")

# qual_scores_before_trim
df3 = pd.read_csv(args.qual_scores_before_trim, delim_whitespace=True, index_col=None)
df3["x"] = df3["Score"] * df3["readsNum"]
sum_reads_num = sum(df3["readsNum"])
phred_avg_before = sum(df3["x"]) / sum_reads_num
phred_avg_before = "{:.2f}".format(phred_avg_before)

# qual_scores_after_trim
df4 = pd.read_csv(args.qual_scores_after_trim, delim_whitespace=True, index_col=None)
df4["x"] = df4["Score"] * df4["readsNum"]
sum_reads_num = sum(df4["readsNum"])
phred_avg_after = sum(df4["x"]) / sum_reads_num
phred_avg_after = "{:.2f}".format(phred_avg_after)

# Preparing output list
output_string = ""
output_list = [
    sample_name,
    reads_before_trim,
    str(GC_content_before),
    str(phred_avg_before),
    reads_after_trim,
    paired_reads_after_trim,
    unpaired_reads_after_trim,
    str(GC_content_after),
    str(phred_avg_after),
    str(coverage),
]

# Creating tab delimited string for qc report generating
for item in output_list:
    output_string += str(item) + "\t"
print(output_string)
