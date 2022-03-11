import sys
import re
import argparse

# take input
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('infile', type=argparse.FileType('r'), nargs='*')

args = parser.parse_args()

args = parser.parse_args()
for f in args.infile:
    print(f)
    # for line in f:
    # process file...
