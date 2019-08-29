import pandas as pd
from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.model_selection import GridSearchCV


class GeneClassifier:

    def __init__(self, input_df, causal_genes, extended_list):
        self.input_df = input_df
        self.causal_genes = causal_genes
        self.extended_list = extended_list

    def return_probability(self, input_df, causal_genes, extended_list):

        x, y = return_x_y(input_df, causal_genes)

        X_train, X_test, y_train, y_test = train_test_split(x, y,
                                                            test_size=0.35,
                                                            random_state=41)

        if sum(y_train) / len(y_train) <= 0.3:
            new_X_train, new_y_train = return_SMOTE(X_train, y_train)

        rf = RandomForestClassifier(random_state=42)

        parameters = best_parameters(rf, new_X_train, new_y_train, PARAM_DIST)

        rf.set_params(bootstrap=True,
                      criterion=parameters['criterion'],
                      max_depth=parameters['max_depth'],
                      max_features=parameters['max_features'],
                      min_samples_leaf=1,
                      min_samples_split=2)


def return_Y(df, db):
    genes = list(set(df['gene_symbol']) & set(db['gene_symbol']))
    y = list(df['gene_symbol'].isin(genes).astype(int))
    print(f'We found {sum(y)} causal genes from this list!')
    return y


def return_gene_symbols(df, db):
    genes = list(set(df['gene_symbol']) & set(db['gene_symbol']))
    return list(df[df['gene_symbol'].isin(genes)]['gene_symbol'])


def return_x_y(gtex_frames, scz_genes):
    y = gtex_frames['gene_symbol'].map(lambda x: 1 if x in scz_genes else 0)
    gtex_frames.drop('gene_id', axis=1, inplace=True)
    x = gtex_frames.fillna(0).set_index('gene_symbol')
    return x, y


def return_SMOTE(X_train, y_train, ratio=0.5):
    sampler = SMOTE(ratio=0.5, random_state=42)
    X_rs, y_rs = sampler.fit_sample(X_train, y_train)
    return X_rs, y_rs


def best_parameters(model, new_X_train, new_y_train, param_dist):
    cv_rf = GridSearchCV(rf, cv=5,
                         param_grid=param_dist,
                         n_jobs=3)
    cv_rf.fit(new_X_train, new_y_train)
    # print('Best Parameters using grid search: \n',
    #   cv_rf.best_params_)
    return cv_rf.best_params_


def error_rate_dict(new_X_train, new_y_train, min_estimators=15, max_estimators=1000):
    error_rate = {}
    for i in range(min_estimators, max_estimators + 1):
        rf.set_params(n_estimators=i)
        rf.fit(new_X_train, new_y_train)
        oob_error = 1 - rf.oob_score_
        error_rate[i] = oob_error
    return error_rate
