import pandas as pd
import os


def get_expression_data(df, gtex_db, gtex_columns):
    """
    This function get median gene expression values for 53 tissues from gtex db
    """
    new_df = pd.DataFrame()
    new_df = new_df.assign(**gtex_columns)
    for index in df.index:
        gene_name = df.loc[index]['gene_symbol']
        try:
            values = gtex_db[gtex_db['Description'] == gene_name].iloc[:, 2:]
            if len(values) == 1:
                new_df = new_df.append(values, ignore_index=True)
            else:
                values = pd.Series(gtex_columns)
                new_df = new_df.append(values, ignore_index=True)
        except:
            print('gtex exception!')
            continue
    return new_df


def transform_to_ranks(df, gtex_columns):
    """
    This function transform median gene expression values to ranks.
    The lowest value will be 0 and the highest expression value - 53.
    """
    import math
    gtex_col = list(gtex_columns.keys())
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


def add_gtex_data(df_with_postgap_data, gtex_columns, gtex_db):
    """
    Combine df with ranks of expression with altered postgap dataframe.
    """
    df_wth_expression_data = get_expression_data(
        df_with_postgap_data, gtex_db, gtex_columns)
    transformed_expression_df = transform_to_ranks(
        df_wth_expression_data, gtex_columns)
    combined_df = pd.concat([df_with_postgap_data.reset_index(
        drop=True), transformed_expression_df], axis=1)
    return combined_df
