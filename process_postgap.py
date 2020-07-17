#!/usr/bin/env python

import sys
import os
import argparse
import glob
import warnings
import pandas as pd
from pathlib import Path
import dask.dataframe as dd
from gprior.var import NEED
from gprior.intoolbox import give_path, merge
from gprior.AdditionalFeatures import add_nsnp  

warnings.filterwarnings("ignore", category=FutureWarning)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='./process_postgap.py',
                            usage='%(prog)s [options]',
                            description="""
                            This script helps to aggregate postgap
                            output for subsequent GPrior usage.
                            """)
    parser.add_argument('-i', '--input', required=True, metavar='',
                        help='Path to file/folder with postgap output')
    parser.add_argument('-o', '--output', required=True, metavar='',
                        help='Path to output')
    
    args = parser.parse_args()
    
    tmp = dd.read_csv(give_path(args.input), sep='\t',
                                             usecols = NEED,
                                             assume_missing=True)
    print('Merging...')
    df = merge(tmp)
    df.columns = ['_'.join(col).strip() for col in df.columns.values]
    df['nSNP'] = add_nsnp(df, tmp)   
    df.to_csv(args.output, sep='\t')
    sys.exit("Done!")
