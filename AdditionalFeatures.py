import pandas as pd
from progress.spinner import Spinner
import math
import os
import sys

def add_reactome_feature(df, db, causal_genes):
    feature = list()
    for gene in df.index:
        if gene in db.Gene1.values:
            n = db.query(f'Gene1 == "{gene}"')['Gene2'] \
             .map(lambda x: 1 if x in causal_genes.gene_symbol.values else 0) \
             .sum()
            feature.append(n)
        else:
            feature.append(0)
    return feature

def get_expression_data(df, GTEX_DB, GTEX_COLUMNS):
    """
    This function get median gene expression values for 53 tissues from gtex db
    """
    new_df = pd.DataFrame().assign(**GTEX_COLUMNS)
    spinner = Spinner('Adding fetures ')
    for gene_name in df.index:
        spinner.next()
        values = GTEX_DB[GTEX_DB['Description'] == gene_name].iloc[:, 2:]
        if len(values) == 1:
            new_df = new_df.append(values, ignore_index=True)
        else:
            values = pd.Series(GTEX_COLUMNS)
            new_df = new_df.append(values, ignore_index=True)
    spinner.finish()
    return new_df


def transform_to_ranks(df, GTEX_COLUMNS):
    """
    This function transform median gene expression values to ranks.
    The lowest value will be 0 and the highest expression value 53.
    """
    gtex_col = list(GTEX_COLUMNS.keys())

    for gene_name in df.index:
        try:
            row_list = [x if math.isnan(
                x) == False else 0 for x in df[gtex_col].loc[gene_name]]
            ranks = {value: rank for rank,
                     value in enumerate(sorted(set(row_list)))}
            ranked = [ranks[i] for i in row_list]
            df.loc[gene_name, gtex_col] = ranked

        except:
            continue
    return df


def add_gtex_feature(df_with_postgap_data, GTEX_COLUMNS, GTEX_DB):
    """
    Combine df with ranks of expression with altered postgap dataframe.
    Add new 53 GTEx features.
    """
    df_wth_expression_data = get_expression_data(
        df_with_postgap_data, GTEX_DB, GTEX_COLUMNS)
    transformed_expression_df = transform_to_ranks(
        df_wth_expression_data, GTEX_COLUMNS)
    transformed_expression_df.index = df_with_postgap_data.index 
    return pd.concat([df_with_postgap_data, transformed_expression_df], axis=1)


def add_gene_similarity_feature(df, db, causal_genes):
    """
    For each gene in df counts how many similar genes are in true_genes.
    Use in order to add several new features:
    1) gtex_similarity
    2) blastp_similarity
    3) atlas_similarity
    4) gene_intaractions
    see: process_input_file function
    """
    feature = list()
    for gene in df.index:
        if gene in db.index:
            gene_list = db.loc[gene, 'associated_genes']
            if type(gene_list) != float:
                gene_list = gene_list.split(",")
                s = sum(
                    list(map(lambda x: x in causal_genes.gene_symbol.values, gene_list)))
                feature.append(s)
            else:
                feature.append(0)
        else:
            feature.append(0)
    return feature
