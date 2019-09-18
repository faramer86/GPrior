import pandas as pd
import os


def get_expression_data(df, GTEX_DB, GTEX_COLUMNS):
    """
    This function get median gene expression values for 53 tissues from gtex db
    """
    new_df = pd.DataFrame()
    new_df = new_df.assign(**GTEX_COLUMNS)
    for index in df.index:
        gene_name = df.loc[index]['gene_symbol']
        try:
            values = GTEX_DB[GTEX_DB['Description'] == gene_name].iloc[:, 2:]
            if len(values) == 1:
                new_df = new_df.append(values, ignore_index=True)
            else:
                values = pd.Series(GTEX_COLUMNS)
                new_df = new_df.append(values, ignore_index=True)
        except:
            print('gtex exception!')
            continue
    return new_df


def transform_to_ranks(df, GTEX_COLUMNS):
    """
    This function transform median gene expression values to ranks.
    The lowest value will be 0 and the highest expression value - 53.
    """
    import math
    gtex_col = list(GTEX_COLUMNS.keys())
    for index in df.index:
        try:
            row_list = [x if math.isnan(
                x) == False else 0 for x in df[gtex_col].loc[index]]
            ranks = {value: rank for rank,
                     value in enumerate(sorted(set(row_list)))}
            ranked = [ranks[i] for i in row_list]
            df.loc[index, gtex_col] = ranked
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
    combined_df = pd.concat([df_with_postgap_data.reset_index(
        drop=True), transformed_expression_df], axis=1)
    print('53 GTEx gene expression features are added')
    return combined_df


def add_gene_similarity_feature(df, db, causal_genes):
    """
    For each gene in df counts how many similar genes are in true_genes.
    Use in order to add two new features: gtex_similarity, blastp_similarity.
    """
    feature = list()
    for gene in df['gene_symbol']:
        if gene in db.index:
            gene_list = db.loc[gene, 'associated_genes']
            if type(gene_list) != float:
                gene_list = gene_list.split(",")
                s = sum(
                    list(map(lambda x: x in causal_genes['gene_symbol'].values, gene_list)))
                feature.append(s)
            else:
                feature.append(0)
        else:
            feature.append(0)
    print('One of the UCSC gene similarity feature is added')
    return feature
