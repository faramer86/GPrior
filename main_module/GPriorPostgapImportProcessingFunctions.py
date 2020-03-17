from main_module.var import *
from progress.spinner import Spinner
import pandas as pd
import os
import glob
from progress.bar import Bar
import main_module.AdditionalFeatures as proc_fts

def drop_columns(df, bar, DO_NOT_NEED):
    """
    Drop unnecessary columns from postgap output df
    """
    bar.next()
    df.drop(DO_NOT_NEED, axis=1, inplace=True)
    return df

def import_postgap_file(input_path, COL_DTYPES):
    """
    Import file of every file in repository.
    """
    if os.path.isfile(input_path):
        with Bar('Importing', max=1) as bar:
            infile = pd.read_csv(input_path, sep='\t', low_memory=False)   # low_memory=False
            if DO_NOT_NEED[0] in list(infile):
                bar.next()
                return drop_columns(infile, DO_NOT_NEED)
            else:
                bar.next()
                return infile
    else:
        try:
            all_files = glob.glob(os.path.join(input_path, "*.tsv"))
        except ZeroDivisionError:
            print('Please, check file extension. Have to be *.tsv')
        with Bar('Importing', max=len(list(all_files))) as bar:
            try:
                return pd.concat(
                    drop_columns(
                    pd.read_csv(filename, sep='\t', low_memory=False),
                    bar, DO_NOT_NEED)
                    for filename in all_files)
            except KeyError:
                print('Importing...')
                return pd.concat(
                    pd.read_csv(filename, sep='\t', low_memory=False)
                    for filename in all_files)

def combine_duplicates(input_df):
    """
    Compute mean values of each feature for every gene name
    """
    nSNP = input_df.groupby(by='gene_symbol').count().iloc[:,0]
    input_df_1 = input_df.groupby(by='gene_symbol').agg(AGG)
    input_df_2 = input_df.groupby(by='gene_symbol').agg(AGG_mean)
    input_df_2.columns = AGG_mean_names
    out_df = pd.concat([input_df_1, input_df_2], axis=1) 
    out_df['nSNP'] = list(nSNP)
    return out_df

def combine_merged(input_df):
    """
    Compute mean values of each feature for every gene name
    """
    out_df = input_df.groupby(by='gene_symbol').mean()
    out_df['nSNP'] = input_df.groupby(by='gene_symbol').agg({'nSNP':'mean'})
    return out_df

def import_merged_files(input_path):
    """
    Aggregate merged files.
    """
    all_files = glob.glob(os.path.join(input_path, "*.tsv"))
    concat = pd.concat([pd.read_csv(file, sep='\t') for file in all_files])
    print(f'There were {len(concat)} gene names.')
    combined = combine_merged(concat) 
    print(f'Now there are {len(combined)} unique gene names.')
    return combined


def process_input_file(postgap_file):
    ml_input = combine_duplicates(postgap_file.set_index('gene_symbol').fillna(0))
    return ml_input
