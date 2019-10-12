from imblearn.over_sampling import SMOTE
from collections import Counter
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.model_selection import KFold, cross_val_score
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from var import *
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


class GeneClassifier:

    def __init__(self, input_df, causal_genes, model, class_weights):
        self.input_df = input_df
        self.causal_genes = causal_genes['gene_symbol'].values
        self.class_weights = class_weights
        self.X_train = 0
        self.X_test = 0
        self.y_train = 0
        self.y_test = 0

    @staticmethod
    def filter_data(df):
        """
        Filter out genes with more then 70% missing features.
        """
        return df.loc[gtex_frames.apply(lambda x: sum(x == 0), axis=1) <= 80, ]

    def return_x_y(self):
        """
        Return X - df with filled name
               y - Series with supervided answers (1 or 0) for each gene
        """
        data = filter_data(self.input_df)
        y = self.input_df['gene_symbol'].map(
            lambda x: 1 if x in self.causal_genes else 0)
        print(f'We found {sum(y)} causal genes from this list!')
        x = self.input_df.fillna(0).set_index('gene_symbol')
        return x, y

    @staticmethod
    def best_parameters(model, X_train, y_train, params):
    cv = GridSearchCV(model, cv=5,
                      param_grid=params,
                      n_jobs=3,
                      iid=False)
    cv.fit(X_train, y_train)
    print('Best Parameters using grid search: \n',
          cv.best_params_)
    return cv.best_params_

    @staticmethod
    def mean_model_prob(plot_df):
        models = ['lr', 'rf', 'dt', 'svc']
        mean_prob = np.array([0.0] * len(plot_df['rf_results']))
        for model in models:
            mean_prob += np.array(plot_dfs[f'{model}_results']['predicted_y'])
        result = pd.DataFrame({'predicted_y': mean_prob / len(models)})
        plot_df['mean_results'] = result
        return plot_df

    @staticmethod
    def PU_learning_result(estimator, x, y, n_estimators=100):

        iP = y[y > 0].index
        iU = y[y <= 0].index

        num_oob = pd.DataFrame(np.zeros(shape=y.shape), index=y.index)
        sum_oob = pd.DataFrame(np.zeros(shape=y.shape), index=y.index)

        for _ in range(n_estimators):
            # Get a bootstrap sample of unlabeled points for this round
            ib = np.random.choice(iU, replace=True, size=3 * len(iP))
            # Find the OOB data points for this round
            i_oob = list(set(iU) - set(ib))
            # Get the training data (ALL positives and the bootstrap
            # sample of unlabeled points) and build the tree
            Xb = x[list(y > 0)].append(x.iloc[ib])
            yb = y[y > 0].append(y.iloc[ib])
            estimator.fit(Xb, yb)
            # Record the OOB scores from this round
            sum_oob.loc[i_oob,
                        0] += estimator.predict_proba(x.iloc[i_oob])[:, 1]
            num_oob.loc[i_oob, 0] += 1
        results = pd.DataFrame(x.index)
        results['Average_prob'] = sum_oob / num_oob
        results = results.sort_values(
            by=['Average_prob'], ascending=False).dropna()
        return results

    def return_probability(self):

        x, y = self.return_x_y()
        print('Splitting data for train/test subsets...')
        self.X_train, self.X_test, self.y_train, self.y_test = \
            train_test_split(x, y, test_size=0.35, random_state=41)

        # if sum(y_train) / len(y_train) <= 0.3:
        #     print('Applying SMOTE...')
        #     print(
        #         f'Classes ration before SMOTE = [1:{sum(y_train)},0:{len(y_train) - sum(y_train)}]')
        #     new_X_train, new_y_train = self.return_SMOTE(X_train, y_train)

        print('Create Classifier and hyperparameters tuning...')

        estimator = set_model_parameters(get_model(self.model))

        parameters = self.best_parameters(
            rf, new_X_train, new_y_train, PARAM_DIST)
        print('Looking for optimal number of trees...')
        n_parameter = self.best_parameters(
            rf, new_X_train, new_y_train, N_ESTIM_PARAM)
        print('Setting parameters to model...')
        rf.set_params(bootstrap=True,
                      criterion=parameters['criterion'],
                      max_depth=parameters['max_depth'],
                      max_features=parameters['max_features'],
                      n_estimators=n_parameter['n_estimators'],
                      min_samples_leaf=1,  # parameters['min_samples_leaf']
                      min_samples_split=2,  # parameters['min_samples_split']
                      class_weight={1: 1, 0: 100})
        rf.fit(new_X_train, new_y_train)
        y_pred = rf.predict(X_test)
        print(confusion_matrix(y_test, y_pred))
        predictions = self.PU_learning_result(rf, x, y)
        return predictions
