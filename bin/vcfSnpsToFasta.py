#!/usr/bin/env python

import sys, re
import argparse
import vcfTools

parser = argparse.ArgumentParser()
parser.add_argument('infile', help='VCF file', type=str)
parser.add_argument('--max_amb_samples', help='maximum number of samples with ambiguous calls for a site to be included', type=int)
parser.add_argument('--min_depth', default=0, help='Replace SNP with "N" if depth is less than minimum (Default: do not check read depth)', type=int)
args = parser.parse_args()

infile = args.infile

max_amb = 1000000
min_depth = args.min_depth
if args.max_amb_samples:
	max_amb = args.max_amb_samples

comment_pattern = re.compile(r"^#")

ref_bases = {}
alt_bases = {}
passed_snp_positions = {}
genome_list = []
amb_pos_counts = {}
sample_translations = {}
line_number = 1

fof_line = infile
sys.stderr.write("Searching " + fof_line)
fof_line = fof_line.rstrip()
header = vcfTools.VcfHeader(fof_line)
caller = header.get_caller()
samples = header.get_samples()
contigs = header.get_contigs()
if line_number == 1:
	for contig in contigs:
		passed_snp_positions[contig] = {}
		amb_pos_counts[contig] = {}
if samples == ['SAMPLE']:
	samples = [fof_line]
	sample_translations[fof_line] = 'SAMPLE'
	sys.stderr.write("No sample name in " + fof_line + ", using file name.\n")
per_file_genome_list = []
for sample in samples:
	genome_list.append(sample)
	per_file_genome_list.append(sample)

with open(fof_line, 'r') as vcf_file:
	for vcf_line in vcf_file:
		if not (re.search(comment_pattern, vcf_line)):
			record = vcfTools.VcfRecord(vcf_line)
			pass_or_fail = record.is_passing(caller)
			for genome in per_file_genome_list:
				translation = genome
				try:
					translation = sample_translations[genome]
				except:
					pass
				genotype = record.get_genotype(index=header.get_sample_index(translation),min_gq=0)
				record_format = dict(zip(record.format.split(':'), record.genotypes[header.get_sample_index(translation)].split(':')))
				variant_type = record.get_variant_type(caller,genotype)
				### print(genome + " " + str(genotype) + " " + str(pass_or_fail) + " " + str(variant_type)) ###
				if pass_or_fail and not variant_type:
					pass
				elif pass_or_fail and variant_type == 'SNP' and int(record_format['DP']) >= min_depth:
					chrom = record.get_chrom()
					pos = int(record.get_pos())

					if not genome in alt_bases:
						alt_bases[genome] = {}
					if not chrom in alt_bases[genome]:
						alt_bases[genome][chrom] = {}
					alt_bases[genome][chrom][pos] = record.get_alt(genotype)
					### print(alt_bases[genome][chrom][pos]) ###
					### print(alt_bases) ###

					if not chrom in ref_bases:
						ref_bases[chrom] = {}
					ref_bases[chrom][pos] = record.get_ref()

					if alt_bases[genome][chrom][pos] != 'N':
						passed_snp_positions[chrom][pos] = True
				else:
					chrom = record.get_chrom()
					pos = int(record.get_pos())
					if not genome in alt_bases:
						alt_bases[genome] = {}
					if not chrom in alt_bases[genome]:
						alt_bases[genome][chrom] = {}
					alt_bases[genome][chrom][pos] = 'N'

					if not pos in amb_pos_counts[chrom]:
						amb_pos_counts[chrom][pos] = 1
					else:
						amb_pos_counts[chrom][pos] += 1
					### print(alt_bases[genome][chrom][pos]) ###
					### print(alt_bases) ###
line_number += 1

### print(alt_bases) ###

sorted_chroms = sorted(passed_snp_positions.keys())

print(">reference")
sequence = ''
for sorted_chrom in sorted_chroms:
	sorted_positions = sorted(passed_snp_positions[sorted_chrom])
	for sorted_position in sorted_positions:
		if sorted_position in amb_pos_counts[sorted_chrom]:
			if amb_pos_counts[sorted_chrom][sorted_position] <= max_amb:
				sequence += ref_bases[sorted_chrom][sorted_position]
		else:
			sequence += ref_bases[sorted_chrom][sorted_position]
for i in range(0,len(sequence),60):
	print(sequence[i:i+60])

num_masked = 0

for ind_genome in genome_list:
	print(">" + ind_genome)
	sequence = ''
	for sorted_chrom in sorted_chroms:
		sorted_positions = sorted(passed_snp_positions[sorted_chrom])
		for sorted_position in sorted_positions:
			if sorted_position in amb_pos_counts[sorted_chrom]:
				if amb_pos_counts[sorted_chrom][sorted_position] <= max_amb:
					try:
						if sorted_position in alt_bases[ind_genome][sorted_chrom]:
							sequence += alt_bases[ind_genome][sorted_chrom][sorted_position]
						else:
							sequence += ref_bases[sorted_chrom][sorted_position]
					except:
						sequence += ref_bases[sorted_chrom][sorted_position]
				else:
					num_masked +=1
			else:
				try:
					if sorted_position in alt_bases[ind_genome][sorted_chrom]:
						sequence += alt_bases[ind_genome][sorted_chrom][sorted_position]
					else:
						sequence += ref_bases[sorted_chrom][sorted_position]
				except:
					sequence += ref_bases[sorted_chrom][sorted_position]
	for i in range(0,len(sequence),60):
		print(sequence[i:i+60])

#sys.stderr.write(str(num_masked) + " sites filtered.\n")
