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
        self.model = model
        self.class_weights = class_weights
        self.X_train = 0
        self.X_test = 0
        self.y_train = 0
        self.y_test = 0

    def return_x_y(self):
        y = self.input_df['gene_symbol'].map(
            lambda x: 1 if x in self.causal_genes else 0)
        print(f'We found {sum(y)} causal genes from this list!')
        x = self.input_df.fillna(0).set_index('gene_symbol')
        return x, y

    def get_model(self):
        if self.model == 'rf':
            print('Random Forest is chosen')
            estimator = RandomForestClassifier(class_weight=self.class_weights)
        elif self.model == 'lr':
            print('Logistic Regression is chosen')
            estimator = LogisticRegression(class_weight=self.class_weights)
        elif self.model == 'xgb':
            print('Extreme Gradient Boosting is chosen')
            estimator = XGBClassifier(
                eval_metric='auc',
                learning_rate=0.1,
                nthread=4,
                silent=True,
                objective='binary:logistic',
                class_weight=self.class_weights)
        print(f'Class weights are: {self.class_weights}')
        return estimator

    def set_model_parameters(self, estimator):
        if self.model == 'rf':
            parameters_rf = best_parameters(
                estimator, self.X_train, self.y_train, PARAM_DIST_RF)
            estimator.set_params(bootstrap=True,
                                 criterion=parameters_rf['criterion'],
                                 max_depth=parameters_rf['max_depth'],
                                 max_features=parameters_rf['max_features'],
                                 min_samples_leaf=1,
                                 min_samples_split=2,
                                 n_estimators=700)
        elif self.model == 'lr':
            parameters_lr = best_parameters(
                lr, self.X_train, self.y_train, PARAM_DIST_LR)
            lr.set_params(C=parameters_lr['C'],
                          penalty=parameters_lr['penalty'])
        elif self.model == 'xgb':
            estimator = XGBClassifier(
                eval_metric='auc',
                learning_rate=0.1,
                nthread=4,
                silent=True,
                objective='binary:logistic',
                class_weight=self.class_weights)
        print('Hyperparameters are tuned..')
        return estimator

    @staticmethod
    def return_SMOTE(X_train, y_train, ratio=0.5):
        print(f'SMOTE ration is {ratio}')
        sampler = SMOTE(ratio=0.5, random_state=42)
        X_rs, y_rs = sampler.fit_sample(X_train, y_train)
        return X_rs, y_rs

    @staticmethod
    def best_parameters(model, new_X_train, new_y_train, PARAM_DIST):
        cv_rf = GridSearchCV(model, cv=5,
                             param_grid=PARAM_DIST,
                             n_jobs=3)
        cv_rf.fit(new_X_train, new_y_train)
        print('Best Parameters using grid search: \n', cv_rf.best_params_)
        return cv_rf.best_params_

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
