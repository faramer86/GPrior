from var import *
from AdditionalFeatures import *
import pandas as pd
import os
import glob


def import_postgap_file(input_path, COL_DTYPES):
    """
    Import file of every file in repository.
    """
    if os.path.isfile(input_path):
        print(f'postgap file "{os.path.basename(input_path)}" is importing..')
        input = drop_columns(pd.read_csv(
            input_path, sep='\t', dtype=COL_DTYPES), DO_NOT_NEED)   # low_memory=False
    else:
        all_files = glob.glob(os.path.join(input_path, "*.txt"))
        li = []
        count = 0
        for filename in all_files:
            if count == 0:
                print(
                    f'postgap file "{os.path.basename(filename)}" from repo is importing..')
                df = drop_columns(pd.read_csv(
                    filename, sep='\t', dtype=COL_DTYPES), DO_NOT_NEED)
                li = df.copy()
                # df.to_csv('tmp.tsv', sep='\t', index=False)
                count += 1
            else:
                print(
                    f'postgap file "{os.path.basename(filename)}" from repo is importing..')
                df = drop_columns(pd.read_csv(
                    filename, sep='\t', dtype=COL_DTYPES), DO_NOT_NEED)    # dtype=COL_DTYPES
                li = pd.concat([li, df.copy()], axis=0,
                               ignore_index=True, sort=False)
                if len(li) >= 10**7:
                    print(
                        'You have reached 10 mil rows of data. \
                        We recommend you to use server or \
                        split your data to chunks and then use this tool!')
                # df.to_csv('tmp.tsv', mode='a', columns=False,
                # index = False, sep = '\t')
                # del df
        # li.append(df)
        print('Combining files...')
        # input = pd.read_csv('tmp.tsv', sep='\t')
        # input = pd.concat(li, axis=0, ignore_index=True, sort=False)
        print('All postgap files from repo are combined')
    return li


def drop_columns(df, DO_NOT_NEED):
    """
    Drop unnecessary columns from postgap output df
    """
    df.drop(DO_NOT_NEED, axis=1, inplace=True)
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
            #gene_name = list(set(input_df.loc[ind]['gene_id']))[0]
            new_df.loc[ind] = mean
            #new_df.loc[ind, 'gene_id'] = gene_name
    new_df = new_df.reset_index()
    new_df.rename(columns={'index': 'gene_symbol'}, inplace=True)
    return new_df


def process_input_file(postgap_file, causal_genes, GTEX_COLUMNS, GTEX_DB):
    """
    Process input postgap file:
    1) Drop unnecessary drop_columns
    2) Fill NA values with zeros
    3) Combine duplicated genes (compute mean values)
    4) Add gtex gene expression ranks for 53 tissues
    5) Add several new features:
        gtex_similarity,
        blastp_similarity,
        gene_interactions
    """
    ML_input = combine_duplicates(fill_NA(postgap_file))
    ML_input = add_gtex_feature(ML_input, GTEX_COLUMNS, GTEX_DB)
    ML_input['gtex_similarity'] = add_gene_similarity_feature(
        ML_input, GTEX_SIMILARITY_DB, causal_genes)
    ML_input['blastp_similarity'] = add_gene_similarity_feature(
        ML_input, BLASTP_SIMILARITY_DB, causal_genes)
    ML_input['gene_interactions'] = add_gene_similarity_feature(
        ML_input, GENE_INTERACTIONS_DB, causal_genes)
    ML_input.drop('gene_id', axis=1, inplace=True)
    return ML_input
