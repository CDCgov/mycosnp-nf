#!/usr/bin/env python3
import argparse
import sys
import re
import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('File_list', metavar='file', nargs='+', type=str)
    
  
    return parser.parse_args(args)
   
def read_data(input_file):
    df=pd.read_csv(input_file)
    return df



def main(args=None):
    args = parse_args(args)
    file_info=[]
    
    clades=['CI','CII','CIII', 'CIV', 'CV', 'CVI']
    #put loop here for files
    for i in args.File_list:
        file_df=read_data(i)
        info=file_df.loc[0, :].values.tolist()
        #list format is 'sample', 'ani', 'target'
        sample=info[0]
        ani=info[1]
        best_match=info[2]
        species='NA'
        clade='NA'
        mash_dist='NA'
        flag='NA'

        for cl in clades:
            pattern=re.compile(rf'\b{re.escape(cl)}\b')
            if (pattern.search(best_match)) and (ani >= 99.7):
                species='Candida_auris'
                clade=cl
                break
            elif (pattern.search(best_match)) and (99.7>ani>95.0):
                species='Candida_auris*'
                clade="UNKNOWN"
                flag="CHECK SAMPLE - highest match to C.auris but low ANI"
                break
            else:
                pass

        if clade=='NA':        
            spec=re.compile(r'\.(.*?)\.')
            match=spec.findall(best_match)
        
            if match:
                if (ani >= 95.0):
                    species=match[0]
                    flag='CHECK_SAMPLE-is not C.auris'
                    
                else:
                    species='UNKNOWN'
                    flag='CHECK SAMPLE - could not determine species'
                    
            else:
                flag='ERROR - could not match any species name in signature files'
        else:
            pass

                
        file_info.append([sample, species, clade, mash_dist, ani, best_match, flag])
    report_df = pd.DataFrame(file_info, columns= ['sample', 'species', 'clade', 'mash_dist', 'ani', 'best_match', 'flag'])
    report_df.to_csv('final_report_mqc.csv', index=False)


if __name__=='__main__':
    sys.exit(main())
