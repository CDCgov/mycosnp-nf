#! /usr/bin/env python3
from lib2to3.pgen2.token import NUMBER
import os
import sys
import pandas as pd
# read input
# Look into collection of read in them
# Open *stats.txt file

# initialize variables
# read length before trimming
read_length_before_trim = 0
# reads before trimming and add to array
reads_before_trim = 0
# reads after trimming and add to array
reads_after_trim = 0
# Paired reads after trimming
paired_reads_after_trim = 0
# unpaired reads after trimming
unpaired_reads_after_trim = 0
# output the numerator (reads count * reads avg length)  for the coverage calculation to be divided by reference length in the module
# coverage_numer = reads_after_trim * #read_avg_length
with open(os.path.join(os.path.dirname(__file__), "/home/tdotrang/qc_reports/ERR2172266.stats.txt"), "r") as f:
    lines = f.readlines()

    reads_before_trim = lines[1].split(" ")[2]
    print(reads_before_trim)

    read_length_before_trim = lines[3].split(" ")[2]
    print(read_length_before_trim)

    reads_after_trim = lines[6].split(" ")[2]
    print(reads_after_trim)

    paired_reads_after_trim = lines[9].split(" ")[5]
    print(paired_reads_after_trim)

    unpaired_reads_after_trim = lines[11].split(" ")[5]
    print(unpaired_reads_after_trim)

coverage_numer = float(read_length_before_trim) * float(reads_before_trim)
print(coverage_numer)

f.close()
# Close *stats.txt file

# Open qa*base_content.txt

df = pd.read_csv('/home/tdotrang/qc_reports/qa.ERR2172266.base_content.txt',
                 delim_whitespace=True, header=None, index_col=None)
print(df)
df_subset = df[df[0] == "GC"]
df_subset[3] = df[1] * df[2]

print(df_subset)
sum_reads = sum(df_subset[2])
sum_reads_GC_content = sum(df_subset[3])
GC_content = (sum_reads_GC_content) / (sum_reads)
print("GC content", GC_content)

# parse only GC lines and calculate average GC content
# Collect the percentages of GC content
# Collect the frequencey of reads with specific GC content
# Multply the percenatage columns by the frequency columns
# Sum the frequencies
# Sum the (percentages*frequncies)
# divide the (sum(percentages)*sum(frequncies)) by total frequencies
# Close *base_content.txt

# Open *base_content.txt
df = pd.read_csv('/home/tdotrang/qc_reports/ERR2172266.base_content.txt',
                 delim_whitespace=True, header=None, index_col=None)
print(df)
df_subset = df[df[0] == "GC"]
df_subset[3] = df[1] * df[2]

print(df_subset)
sum_reads = sum(df_subset[2])
sum_reads_GC_content = sum(df_subset[3])
GC_content = (sum_reads_GC_content) / (sum_reads)
print("GC content", GC_content)
# parse only GC lines and calculate average GC content
# Collect the percentages of GC content
# Collect the frequencey of reads with specific GC content
# Multply the percenatage columns by the frequency columns
# Sum the frequencies
# Sum the (percentages*frequncies)
# divide the (sum(percentages)*sum(frequncies)) by total frequencies
# Close qa*base_content.txt

# THE qa.* FILES ARE THE AFTER DATA!!!!!!!!
# Open *for_qual_histograms.txt

f = pd.read_csv('/home/tdotrang/qc_reports/qa.ERR2172266.for_qual_histogram.txt',
                delim_whitespace=True, index_col=None)

f['x'] = f['Score']*f['readsNum']

sum_reads_num = sum(f['readsNum'])
print("sum of reads: ", sum_reads_num)

phred_avg = sum(f['x'])/sum_reads_num
print("phred quality score: \n", f)
print("phred average: ", phred_avg)

# no pandas methods
# with open(os.path.join(os.path.dirname(__file__), "/home/tdotrang/qc_reports/ERR2172266.for_qual_histogram.txt"), "r") as f:
#     # Collect the PHRED scores
#     phred_score_array = []
#     next(f)
#     for line in f:
#         phred_score_array.append(line.strip().split('\t')[0:2])
#         for i in phred_score_array:
#             phred_freq = [0] * [1]
#             i.append(phred_freq)
#     print(phred_score_array)


# Collect the frequencey of reads with specific PHRED content
# Check if we want read number or read bases?
# Multply the PHRED columns by the frequency columns
# Sum the frequencies
# Sum the (PHRED*frequncies)
# divide the (sum(PHRED)*sum(frequncies)) by total frequencies
# Close *for_qual_histograms.txt

# Open qa*for_qual_histograms.txt
f = pd.read_csv('/home/tdotrang/qc_reports/ERR2172266.for_qual_histogram.txt',
                delim_whitespace=True, index_col=None)

f['x'] = f['Score']*f['readsNum']

sum_reads_num = sum(f['readsNum'])
print("AFTER sum of reads: ", sum_reads_num)

phred_avg = sum(f['x'])/sum_reads_num
print("phred quality score: \n", f)
print("AFTER phred average: ", phred_avg)
# Collect the PHRED scores
# Collect the frequencey of reads with specific PHRED content
# Multply the PHRED columns by the frequency columns
# Sum the frequencies
# Sum the (PHRED*frequncies)
# divide the (sum(PHRED)*sum(frequncies)) by total frequencies
# Close qa*for_qual_histograms.txt
