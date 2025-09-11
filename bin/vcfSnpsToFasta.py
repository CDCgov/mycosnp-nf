#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict
import vcfTools


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='VCF file', type=str)
    parser.add_argument('--max_amb_samples', help='Maximum number of samples with ambiguous calls for a site to be included', type=int, default=1000000)
    parser.add_argument('--min_depth', default=0, help='Replace SNP with "N" if depth is less than minimum (Default: do not check read depth)', type=int)
    return parser.parse_args()


def initialize_data_structures(contigs, samples):
    passed_snp_positions = {contig: {} for contig in contigs}
    amb_pos_counts = {contig: {} for contig in contigs}
    alt_bases = defaultdict(lambda: defaultdict(dict))
    ref_bases = defaultdict(dict)
    genome_list = samples[:]
    return passed_snp_positions, amb_pos_counts, alt_bases, ref_bases, genome_list


def process_vcf_record(record, header, caller, genome, min_depth, alt_bases, ref_bases, passed_snp_positions, amb_pos_counts):
    translation = genome
    genotype = record.get_genotype(index=header.get_sample_index(translation), min_gq=0)
    record_format = dict(zip(record.format.split(':'), record.genotypes[header.get_sample_index(translation)].split(':')))
    variant_type = record.get_variant_type(caller, genotype)
    chrom = record.get_chrom()
    pos = int(record.get_pos())

    if record.is_passing(caller) and variant_type == 'SNP' and int(record_format.get('DP', 0)) >= min_depth:
        alt_bases[genome][chrom][pos] = record.get_alt(genotype)
        ref_bases[chrom][pos] = record.get_ref()
        if alt_bases[genome][chrom][pos] != 'N':
            passed_snp_positions[chrom][pos] = True
    else:
        alt_bases[genome][chrom][pos] = 'N'
        amb_pos_counts[chrom][pos] = amb_pos_counts[chrom].get(pos, 0) + 1


def generate_fasta(genome_list, sorted_chroms, passed_snp_positions, amb_pos_counts, ref_bases, alt_bases, max_amb):
    def build_sequence(genome):
        sequence = ''
        for chrom in sorted_chroms:
            sorted_positions = sorted(passed_snp_positions[chrom])
            for pos in sorted_positions:
                if pos in amb_pos_counts[chrom] and amb_pos_counts[chrom][pos] > max_amb:
                    sequence += 'N'
                else:
                    sequence += alt_bases[genome][chrom].get(pos, ref_bases[chrom][pos])
        return sequence

    print(">reference")
    ref_sequence = build_sequence('reference')
    for i in range(0, len(ref_sequence), 60):
        print(ref_sequence[i:i+60])

    for genome in genome_list:
        print(f">{genome}")
        genome_sequence = build_sequence(genome)
        for i in range(0, len(genome_sequence), 60):
            print(genome_sequence[i:i+60])


def main():
    args = parse_arguments()
    infile = args.infile
    max_amb = args.max_amb_samples
    min_depth = args.min_depth

    sys.stderr.write(f"Searching {infile}\n")
    header = vcfTools.VcfHeader(infile)
    caller = header.get_caller()
    samples = header.get_samples()
    contigs = header.get_contigs()

    if samples == ['SAMPLE']:
        samples = [infile]
        sys.stderr.write(f"No sample name in {infile}, using file name.\n")

    passed_snp_positions, amb_pos_counts, alt_bases, ref_bases, genome_list = initialize_data_structures(contigs, samples)

    with open(infile, 'r') as vcf_file:
        for vcf_line in vcf_file:
            if not vcf_line.startswith('#'):
                record = vcfTools.VcfRecord(vcf_line)
                for genome in genome_list:
                    process_vcf_record(record, header, caller, genome, min_depth, alt_bases, ref_bases, passed_snp_positions, amb_pos_counts)

    sorted_chroms = sorted(passed_snp_positions.keys())
    generate_fasta(genome_list, sorted_chroms, passed_snp_positions, amb_pos_counts, ref_bases, alt_bases, max_amb)


if __name__ == "__main__":
    main()
