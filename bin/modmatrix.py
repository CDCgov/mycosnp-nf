#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("InputFile", type=str)
    parser.add_argument("OutputFile", type=str)
    return parser.parse_args(args)
    

def read_data(input_file):
    anim=pd.read_csv(input_file)
    return anim

def reformat(anim):
    x=list(anim.columns)
    anim.index = x
    anim = anim.mul(100.000000) 
    
    out_df=pd.DataFrame(columns=['source', 'value', 'target'])
    for index, row in anim.iterrows():
        for column, value in row.items():
            if index != column:
                out_df.loc[len(out_df.index)]= [index, value, column]

    out_df = out_df.sort_values(by='source').reset_index(drop=True)
    return out_df
    #return anim

def removeDups(anim):
    unique_comparisons = set()

    for index, row in anim.iterrows():
        comparison=tuple(sorted([row['source'], row['target']]))
        unique_comparisons.add(comparison)
    
    out_anim = pd.DataFrame(columns =['source', 'value', 'target'])
    for comparison in unique_comparisons:
        source, target = comparison
        rows = anim[(anim['source'] == source) & (anim['target'] == target)]
        out_anim = pd.concat([out_anim, rows])

    return out_anim
  
###block


def main(args=None):
    args = parse_args(args)
    ani_m = read_data(args.InputFile)
    ani_ref = reformat(ani_m)
    ani_final=removeDups(ani_ref)

    ani_final.to_csv(args.OutputFile, sep="\t", index=False)

if __name__=='__main__':
    sys.exit(main())