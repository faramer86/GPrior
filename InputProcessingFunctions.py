from var import *
from AdditionalFeatures import *
from progress.spinner import Spinner
import pandas as pd
import os
import glob
from progress.bar import Bar

def number_of_intersections(list_1, list_2):
    return len(set(list_1) & set(list_2))


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
                    pd.read_csv(filename, sep='\t', low_memory=True),
                    bar, DO_NOT_NEED)
                    for filename in all_files)
            except KeyError:
                print('Importing...')
                return pd.concat(
                    pd.read_csv(filename, sep='\t', low_memory=True)
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
    out_df['nSNP'] = input_df.groupby(by='gene_symbol').agg({'nSNP':'max'})
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

def add_features(df, causal_genes):
    df['Reactome'] = add_reactome_feature(
        df, REACTOME_DB, causal_genes)
    df = add_gtex_feature(df, GTEX_COLUMNS, GTEX_DB)
    df['gtex_similarity'] = add_gene_similarity_feature(
        df, GTEX_SIMILARITY_DB, causal_genes)
    df['blastp_similarity'] = add_gene_similarity_feature(
        df, BLASTP_SIMILARITY_DB, causal_genes)
    df['atlas_similarity'] = add_gene_similarity_feature(
        df, ATLAS_SIMILARITY_DB, causal_genes)
    df['gene_interactions'] = add_gene_similarity_feature(
        df, GENE_INTERACTIONS_DB, causal_genes)
    return df


def process_input_file(postgap_file, causal_genes):
    """
    Process input postgap file:
    1) Drop unnecessary drop_columns
    2) Fill NA values with zeros
    3) Combine duplicated genes (compute mean values)
    4) Add gtex gene expression ranks for 53 tissues
    5) Add several new features:
        gtex_similarity,
        atlas_similarity,
        blastp_similarity,
        gene_interactions
    """
    ml_input = combine_duplicates(postgap_file.set_index('gene_symbol').fillna(0))
    print(f'Total number of genes: {len(ml_input.index)}')
    intersect_n = len(set(ml_input.index) & set(causal_genes.gene_symbol))
    print(f'Number of causal genes found: {intersect_n}')
    return add_features(ml_input, causal_genes)
