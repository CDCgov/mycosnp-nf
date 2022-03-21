from genericpath import sameopenfile
from importlib.resources import path
import re
import os
import sys
import argparse
import pandas as pd

# take input
# Separate function to parse args
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample")
parser.add_argument("--stats", type=argparse.FileType("r"))
parser.add_argument("--base_content_before_trim", type=argparse.FileType("r"))
parser.add_argument("--base_content_after_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_before_trim", type=argparse.FileType("r"))
parser.add_argument("--qual_scores_after_trim", type=argparse.FileType("r"))
parser.add_argument("--reference", type=argparse.FileType("r"))

args = parser.parse_args()

# Separate function to produce list or dataframe for metrics, need to read in reference

sample_name = args.sample

# Print
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
f1 = pd.read_csv(args.qual_scores_before_trim, delim_whitespace=True, index_col=None)

f1["x"] = f1["Score"] * f1["readsNum"]

sum_reads_num = sum(f1["readsNum"])

phred_avg_before = sum(f1["x"]) / sum_reads_num
phred_avg_before = "{:.2f}".format(phred_avg_before)


# qual_scores_after_trim
f2 = pd.read_csv(args.qual_scores_after_trim, delim_whitespace=True, index_col=None)

f2["x"] = f2["Score"] * f2["readsNum"]

sum_reads_num = sum(f2["readsNum"])

phred_avg_after = sum(f2["x"]) / sum_reads_num
phred_avg_after = "{:.2f}".format(phred_avg_after)


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

# d = {'Sample Name': sample_name, 'Reads Before Trimming': reads_before_trim , 'GC Before Trimming': GC_content_before, 'Average Phred Before Trimming': phred_avg_before, 'Reads After Trimming': reads_after_trim, 'Paired Reads After Trimming': paired_reads_after_trim, 'Unpaired Reads After Trimming': unpaired_reads_after_trim, 'GC After Trimming': GC_content_after,'Average Phred After Trimming': phred_avg_after,'Coverage After Trimming': coverage}
# df = pd.DataFrame(d)
# df.to_csv("qcreport.tsv", sep="\t")
# print(df)
