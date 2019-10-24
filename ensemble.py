# from sklearn.base import BaseEstimator
# from sklearn.base import ClassifierMixin
# import operator
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.datasets import make_blobs
from progress.bar import IncrementalBar
import pandas as pd
from var import *
import numpy as np
import random
import math


class BootstrapSampleTuner():

    def __init__(self, estimator, params):
        self.estimator = estimator
        self.params = params
    
    def get_tuned_clf(self, Xb, yb):
        cv = GridSearchCV(self.estimator, cv=5, param_grid=self.params, iid=False)
        cv.fit(Xb, yb)
        self.best_params = cv.best_params_
        self.estimator.set_params(**cv.best_params_)
        return self.estimator
    
    def print_best_params(self):
        print(self.best_params)
        pass


class PUBaggingClassifier():

    def __init__(self, clf, params, n_estimators=100, under_sample_coef=0.8):
        self.clf = clf
        self.params = params
        self.n_estimators = n_estimators
        self.under_sample_coef = under_sample_coef

    def run_bagging(self, X, y):

        gene_names = X.index.copy()
        self.gene_names = gene_names

        iP = y[y > 0].index
        self.iP = iP
        iU = y[y <= 0].index
        
        num_oob = pd.DataFrame(np.zeros(shape = y.shape), index = y.index)
        sum_oob = pd.DataFrame(np.zeros(shape = y.shape), index = y.index)
        
        estimator_name = list(self.clf.keys())[0]

        bar = IncrementalBar(f'{estimator_name}', max=self.n_estimators)

        for _ in range(self.n_estimators):
                    
            ib = np.random.choice(iU, replace = True, size = math.ceil(len(iP)*self.under_sample_coef))
            i_oob = list(set(iU) - set(ib))
        
            Xb = X[list(y>0)].append(X.iloc[ib])
            yb = y[y > 0].append(y.iloc[ib])
                    
            estimator_base = list(self.clf.values())[0]            
            estimator = BootstrapSampleTuner(estimator_base, self.params).get_tuned_clf(Xb, yb)            
            estimator.fit(Xb, yb)
            
            sum_oob.loc[i_oob, 0] += estimator.predict_proba(X.iloc[i_oob])[:,1]
            num_oob.loc[i_oob, 0] += 1
            bar.next()

        clf_df = pd.DataFrame({'gene_symbol':gene_names})
        clf_df['Average_prob'] = 100 * sum_oob / num_oob

        self.clf_df = clf_df 
        bar.finish()
        return clf_df
    
    @staticmethod
    def add_y_answers(self, extended_list):
        return self.clf_df['gene_symbol'].map(lambda gene: (gene in extended_list) and \
                                     (gene not in self.gene_names[self.iP].values))
    
    @staticmethod        
    def get_threshold(self, threshold_range, qc_store):
        opt = [threshold_range[i] for i, j in enumerate(qc_store) if j == max(qc_store)][0]
        return opt
    
    @staticmethod
    def calculate_qc_coef(self, threshold):
        self.clf_df['pred_y'] = self.clf_df['Average_prob'] >= threshold
        if sum(self.clf_df['pred_y']) == 0:
            return 0
        Pr = sum(self.clf_df['pred_y']) / len(self.clf_df['pred_y'])
        TP = sum(self.clf_df['y'] & self.clf_df['pred_y']) # TP
        FN = sum(self.clf_df['y']) - TP # FN
        recall = TP/(TP+FN)
        coef = recall**2/Pr
        return coef

    def perform_qc(self, extended_list):
        self.clf_df['y'] = self.add_y_answers(extended_list)
        qc_range = [calculate_qc_coef(self.clf_df, i) for i in THRESHOLD_VALUES]
        optimal_thr = self.get_threshold(threshold_range, qc_range)
        self.clf_df['Optimal_threshold'] = optimal_thr
        self.clf_df['QC_coef'] = self.calculate_qc_coef(optimal)
        return self.clf_df




class EnsembleClassifier():

    def __init__(self, estimators):



