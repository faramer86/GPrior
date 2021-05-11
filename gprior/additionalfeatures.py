import pandas as pd
from gprior.var import *
from progress.spinner import Spinner
import dask.dataframe as dd


def add_nsnp(df, tmp):
    return list(dd.concat([tmp]) \
                .fillna(0) \
                .groupby('gene_symbol') \
                .count() \
                .score \
                .compute())


def add_reactome_feature(df, db, causal_genes):
    feature = list()
    for gene in df.gene_symbol:
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
    for gene_name in df.gene_symbol:
        try:
            spinner.next()
            values = GTEX_DB.query(f'Description == "{gene_name}"').iloc[:, 2:]
            if len(values) != 1:
                values = pd.Series(GTEX_COLUMNS)
            new_df = new_df.append(values, ignore_index=True)
        except:
            continue
    spinner.finish()
    print(new_df)
    return new_df


def add_gtex_feature(df_with_postgap_data, GTEX_COLUMNS, GTEX_DB):
    """
    Combine df with ranks of expression with altered postgap dataframe.
    Add new 53 GTEx features.
    """
    df_wth_expression_data = get_expression_data(
        df_with_postgap_data, GTEX_DB, GTEX_COLUMNS)
    df_wth_expression_data.index = df_with_postgap_data.index
    return pd.concat([df_with_postgap_data, df_wth_expression_data], axis=1)


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
    for gene in df.gene_symbol:
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


def add_features(df, causal_genes):
    df = add_gtex_feature(df, GTEX_COLUMNS, GTEX_DB)
    df['Reactome'] = add_reactome_feature(
        df, REACTOME_DB, causal_genes)
    df['gtex_similarity'] = add_gene_similarity_feature(
        df, GTEX_SIMILARITY_DB, causal_genes)
    df['blastp_similarity'] = add_gene_similarity_feature(
        df, BLASTP_SIMILARITY_DB, causal_genes)
    df['atlas_similarity'] = add_gene_similarity_feature(
        df, ATLAS_SIMILARITY_DB, causal_genes)
    df['gene_interactions'] = add_gene_similarity_feature(
        df, GENE_INTERACTIONS_DB, causal_genes)
    return df
