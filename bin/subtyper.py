#!/usr/bin/env python3

"""
subtyper.py
Authors:
    Jared Johnson, jared.johnson@doh.wa.gov
    Zack Mudge, ZMudge@cdc.gov
"""

import argparse
import csv
import os
import subprocess
import sys

version = "v2.0.0"

def find_highest_ani_sketch(csv_file):
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)[1:]  # Headers line: excluding first column
        values = next(reader)[1:]   # ANI values line: excluding first column

        # Convert ANI values to floats for comparison
        values = list(map(float, values))

        # Find the index of the maximum ANI
        max_index = values.index(max(values))
        highest_sketch = headers[max_index]

        return highest_sketch, values[max_index]

def main():
    parser = argparse.ArgumentParser(description='Subtyper script for the pre-MycoSNP workflow')

    parser.add_argument('sample', help='Sample name')
    parser.add_argument('taxon', help='Taxon name (use underscores for spaces)')
    parser.add_argument('db', help='Path to the subtype database')
    parser.add_argument('seq', help='Path to the assembly or reads')
    parser.add_argument('output_csv', help='Path for the output CSV file')
    parser.add_argument('-v', '--version', action='version', version=version, help='Show version and exit')
    args = parser.parse_args()

    sample = args.sample
    taxon = args.taxon.replace("_", " ")
    db = args.db.rstrip('/')
    seq = args.seq
    output_csv_path = args.output_csv.rstrip('/')

    #----- SUBTYPER -----#
    # Select sourmash signature file (if listed)
    sourmash_taxa_csv = os.path.join(db, "sourmash_taxa.csv")
    if os.path.isfile(sourmash_taxa_csv):
        with open(sourmash_taxa_csv, 'r') as f:
            signature_filepath = None
            for line in f:
                parts = line.strip().split(',')
                if parts[0] == taxon:
                    signature_filepath = parts[1]
                    break
    else:
        print(f"Error: {sourmash_taxa_csv} does not exist.")
        sys.exit(1)

    # Run sourmash if signatures exist
    if signature_filepath and os.path.isfile(os.path.join(db, signature_filepath)):
        # Run sourmash sketch on seq
        sketch_file = f"{sample}.k31.sig"
        subprocess.run(["sourmash", "sketch", "dna", "-p", "k=31", "-o", sketch_file, seq], check=True)

        # Compare sketch of seq to the sketches in the signature file
        ani_comp_csv = f"{sample}.k31.ani_comp.csv"
        subprocess.run(["sourmash", "compare", "-k", "31", "--containment", "--ani", "--csv", ani_comp_csv, sketch_file, os.path.join(db, signature_filepath)], check=True)

        # Find the sketch with the highest ANI
        highest_ani_sketch, ani = find_highest_ani_sketch(ani_comp_csv)
        # Convert ANI to percentage and round to 4 decimal places
        ani_percentage = round(ani * 100, 4)

        # Write the results to the output CSV
        with open(output_csv_path, 'w') as output_csv:
            output_csv.write("sample,subtype_closest_match,est_ANI\n")
            output_csv.write(f"{sample},{highest_ani_sketch},{ani_percentage}\n")

    else:
        print(f"Signature file for {taxon} not found.")

if __name__ == '__main__':
    main()