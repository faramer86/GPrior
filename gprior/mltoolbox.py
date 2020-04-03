
import pandas as pd
from sklearn.cluster import FeatureAgglomeration
from sklearn.preprocessing import RobustScaler
import sys
    
def filter_data(df):
    """
    Filter out genes with more then 85% missing features.
    """
    bool_index = df.apply(lambda x: sum(x == 0), axis=1) <= 85
    return df.loc[bool_index, ]

def process_x(x, n_clusters):
    """
    Normalize and remove multicollineaarity usind FA
    """
    gene_names = x.index
    scaler = RobustScaler()
    try:
        norm_x = scaler.fit_transform(x)
    except ValueError:
        sys.exit('(INPUT) Please, check the requirements for input file!')
    if n_clusters:
        n_clf = n_clusters
    elif len(list(x)) >= 1:
        n_cls = len(list(x))
    else:
        sys.exit('There are no features!')
    return pd.DataFrame(
           FeatureAgglomeration(n_clusters=n_clf).fit_transform(norm_x), 
           index=gene_names)
    

def prepare_df(df):
    """
    Prepare dataframe for RobustScaling and/or
    for Feature Agglomeration
    """
    return df.fillna(0).set_index('gene_symbol')

def prepare_y(df, causal_genes):
    """
    Make Y from provided dataframe and list of causal genes
    We need it for subsequent ML step
    """
    return df.gene_symbol.map(lambda x: 1 if x in causal_genes.gene_symbol.values else 0)

def return_x_y(df, causal_genes, k_clusters):
    """
    Return X - df with filled name
    y - Series with supervided answers (1 or 0) for each gene
    """
    #filtered_df = filter_data(df)
    y = prepare_y(df, causal_genes)
    if sum(y) <= 5:
        raise AssertionError("Numer of causal genes found <= 5. Add more causal genes.")
    x = process_x(prepare_df(df), k_clusters)
    return x, y    
    
