#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("InputFile", type=str)
    parser.add_argument("Sample", type=str)

    return parser.parse_args(args)
   


def filter_sample(comp_file, samplename):
    df = pd.read_csv(comp_file, sep='\t')
    mask = (df['source'].str.contains(samplename, case=False) | df['target'].str.contains(samplename, case=False))
    filt_df=df[mask]

    mask = filt_df['target'].str.contains(samplename, case=False)
    filt_df.loc[mask, ['source','target']] = filt_df.loc[mask,['target', 'source']].values
    filt_df=filt_df.sort_values(by='value', ascending=False)
    return filt_df


def main(args=None):
    args = parse_args(args)
    filt_df = filter_sample(args.InputFile, args.Sample)
    outfile=args.Sample+"_results.csv"
    filt_df.to_csv(outfile, index=False, header=False)
    highest=filt_df.loc[filt_df['value'].idxmax()].values.tolist()
    

if __name__=='__main__':
    sys.exit(main())

