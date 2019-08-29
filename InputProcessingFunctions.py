from var import *
from GtexExpression import *
import pandas as pd
import os
import glob


def drop_columns(df):
    """
    Drop unnecessary columns from postgap output df
    """
    df.drop(do_not_need, axis=1, inplace=True)
    return df


def fill_NA(df):
    """
    Fill NA values for postgap output df.
    It is necessary for subsequent analysis.
    """
    df = df.fillna(0).copy()
    return df


def combine_duplicates(input_df):
    """
    Compute mean values of each feature for every gene name
    """
    input_df = input_df.set_index('gene_symbol')
    new_df = pd.DataFrame(columns=list(input_df))
    for ind in set(input_df.index):
        if input_df.loc[ind].ndim == 1:
            new_df.loc[ind] = input_df.loc[ind]
        else:
            mean = input_df.loc[ind].mean()
            gene_name = list(set(input_df.loc[ind]['gene_id']))[0]
            new_df.loc[ind] = mean
            new_df.loc[ind, 'gene_id'] = gene_name
    new_df = new_df.reset_index()
    new_df.rename(columns={'index': 'gene_symbol'}, inplace=True)
    return new_df


def process_input_file(postgap_file, gtex_columns, gtex_db):
    """
    Process input postgap file:
    1) Drop unnecessary drop_columns
    2) Fill NA values with zeros
    3) Combine duplicated genes (compute mean values)
    4) Add gtex gene expression ranks for 53 tissues
    """
    ML_input = add_gtex_data(combine_duplicates(
        fill_NA(postgap_file)), gtex_columns, gtex_db)
    return ML_input


def import_postgap_file(input_path):
    """
    Depends on path import file of every file in repository
    """
    if os.path.isfile(input_path):
        input = drop_columns(pd.read_csv(
            input_path, sep='\t', dtype=col_dtypes))
        print(f'postgap file "{os.path.basename(input_path)}" is imported')
    else:
        all_files = glob.glob(os.path.join(input_path, "*.txt"))
        li = []
        for filename in all_files:
            df = drop_columns(pd.read_csv(
                filename, sep='\t', dtype=col_dtypes))
            li.append(df)
            print(
                f'postgap file "{os.path.basename(filename)}" from repo is imported')
        print('Combining files...')
        input = pd.concat(li, axis=0, ignore_index=True, sort=False)
        print('All postgap files from repo are combined')
    return input
