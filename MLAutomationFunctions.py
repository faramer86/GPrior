
import pandas as pd
from sklearn.cluster import FeatureAgglomeration
from sklearn.preprocessing import RobustScaler

    
def filter_data(df):
    """
    Filter out genes with more then 80% missing features.
    """
    bool_index = df.apply(lambda x: sum(x == 0), axis=1) <= 80
    return df.loc[bool_index, ]

def process_x(x):
    """
    Normalize and remove feature correlation usind FA
    """
    gene_names = x.copy().index
    scaler = RobustScaler()
    print('Normalizing..')
    norm_x = scaler.fit_transform(x)
    print('Feature aggregation..')
    new_x = pd.DataFrame(FeatureAgglomeration(n_clusters=50).fit_transform(norm_x), index=gene_names)
    return new_x

def return_x_y(df, causal_genes):
    """
    Return X - df with filled name
    y - Series with supervided answers (1 or 0) for each gene
    """
    data = filter_data(df)
    y = df['gene_symbol'].map(
    lambda x: 1 if x in causal_genes['gene_symbol'].values else 0)
    print(f'We found {sum(y)} genes in list of causal genes!')
    if sum(y) <= 3:
        raise AssertionError("Numer of causal genes found <= 3. Add more causal genes.")
    x = process_x(df.fillna(0).set_index('gene_symbol'))
    return x, y    
    
