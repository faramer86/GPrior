import dask.dataframe as dd
import pandas as pd
import os
import sys
import glob

def give_path(arg):
        if os.path.isfile(arg):
            return arg
        else:
            return os.path.join(arg, "*.tsv")

def merge(df):
    return dd.concat([df]) \
             .fillna(0) \
             .groupby('gene_symbol') \
             .agg(['mean', 'max']) \
             .compute()

def find_genes(list_1, list_2):
    return len(set(list_1) & set(list_2))

def give_n(arg, default=15):
    if arg:
        print(f'Number of bootstraps: {arg}')
        return arg
    else:
        print(f'Number of bootstraps: {default}')
        return default

def give_k(arg, input_file):
    if arg:
        n_features = len(list(input_file)) - 1
        print(f'Number of clusters used: {arg}/{n_features}')
        return arg
    else:
        n_features = len(list(input_file)) - 1
        print(f'Number of clusters used: {n_features}')
        return n_features

def give_s(arg, default=1):
    if arg:
        print(f'Sampling coefficient: {arg}\n')
        return arg
    else:
        print(f'Sampling coefficient: {default}\n')
        return default


def import_file(input_path):
    """
    Import files (--input, -ts, -ass)
    """
    filename = os.path.basename(input_path)

    assert '.tsv' in filename, 'Plese, use .tsv extansion for all of the files!'
    
    if os.path.isfile(input_path):
        print(f'Importing {filename}')
        infile = pd.read_csv(input_path, sep='\t', low_memory=False)
        if 'gene_symbol' not in list(infile):
            sys.exit(f'({filename}) Please specify column with genes in all of the files as "gene_symbol"!')
        return infile
    else:
        sys.exit(f'There is no file with such name ({filename}) or it is not a file!')

