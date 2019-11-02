from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.datasets import make_blobs
from progress.bar import IncrementalBar, Bar
from sklearn.base import clone
from itertools import combinations
import pandas as pd
from var import *
import numpy as np
import random
import math

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


class PUBaggingClassifier():

    def __init__(self, clf_name, clf_base, params, n_bootstrap=3001, under_sample_coef=0.7):
        self.clf_name = clf_name
        self.clf_base = clf_base
        self.params = params
        self.n_bootstrap = n_bootstrap
        self.under_sample_coef = under_sample_coef

    def run_bagging(self, X, y):

        self.gene_names = X.index.copy() 
        self.iP = y[y > 0].index
        self.iU = y[y <= 0].index
        
        num_oob = pd.DataFrame(np.zeros(shape = y.shape), index = y.index)
        sum_oob = pd.DataFrame(np.zeros(shape = y.shape), index = y.index)
        
        print(f'{self.clf_name}', end='\t')
        bar = IncrementalBar(max=self.n_bootstrap)
        print()

        for _ in range(self.n_bootstrap):
                    
            ib = np.random.choice(self.iU, replace = True, size = math.ceil(len(self.iP)*self.under_sample_coef))
            i_oob = list(set(self.iU) - set(ib))
        
            Xb = X[list(y>0)].append(X.iloc[ib])
            yb = y[y > 0].append(y.iloc[ib])

            clone_base = clone(self.clf_base)        
            estimator = clone(self.clf_base)
            
            estimator.fit(Xb, yb)
            sum_oob.loc[i_oob, 0] += estimator.predict_proba(X.iloc[i_oob])[:,1]
            num_oob.loc[i_oob, 0] += 1
            bar.next()
                
        clf_df = pd.DataFrame({'gene_symbol': self.gene_names})
        clf_df['Average_prob'] = 100 * sum_oob / num_oob
        bar.finish()
        self.clf_df = clf_df 
        
    
    def add_y_answers(self, extended_list): #['gene_symbol'].values
        return self.clf_df['gene_symbol'].map(lambda gene: (gene in extended_list) & \
                                     (gene not in self.gene_names[self.iP].values))
        
    def get_threshold(self, threshold_range, qc_store):
        return [threshold_range[i] for i, j in enumerate(qc_store) if j == max(qc_store)][0]
        
    def calculate_qc_coef(self, threshold):
        self.clf_df['pred_y'] = self.clf_df['Average_prob'] >= threshold
        if self.clf_df['pred_y'].sum() == 0:
            return 0
        
        Pr = self.clf_df['pred_y'].sum() / len(self.clf_df['pred_y']) 
        
        TP = sum(self.clf_df['y'] & self.clf_df['pred_y']) # TP
        FN = self.clf_df['y'].sum() - TP # FN     
        recall = TP/(TP+FN)        
        coef = recall**2/Pr
        return coef

    def get_df_with_qc(self, extended_list):
        threshold_range = [0.01, 0.1] + [i for i in range(1, 98)]
        self.clf_df['y'] = self.add_y_answers(extended_list)
        qc_range = [self.calculate_qc_coef(i) for i in threshold_range]
        optimal_thr = self.get_threshold(threshold_range, qc_range)
        qc_coef = self.calculate_qc_coef(optimal_thr)
        self.clf_df['Optimal_threshold'] = optimal_thr
        self.clf_df['QC_coef'] = qc_coef
        print(f'QC coefficient: {qc_coef}')
        print(f'Optimal threshold: {optimal_thr}\n')
        return self.clf_df




class EnsembleClassifier(PUBaggingClassifier):

    def __init__(self, X, y, dict_of_estimators, extended_list):
        self.dict_of_estimators = dict_of_estimators
        self.extended_list = set(extended_list['gene_symbol'].values)
        self.X = X
        self.y = y
        self.gene_names = X.index.copy() 
        self.iP = y[y > 0].index.copy()
        self.probas = list()
        self.weights = list()
        self.thresholds = list()

    def wmean_qc_coef(self, wmean, threshold): # ['gene_symbol'].values
        pred_y = np.array(wmean >= threshold)
        new_y = np.array(list(map(lambda gene: (gene in self.extended_list) & \
                        (gene not in self.gene_names[self.iP]), self.gene_names)))
        if np.sum(pred_y) == 0:
            return 0
        Pr = np.sum(pred_y) / pred_y.size
        TP = np.sum(new_y & pred_y) # TP
        FN = np.sum(new_y) - TP # FN
        recall = TP/(TP+FN)
        coef = recall**2/Pr
        return coef
    
    def perform_qc(self, wmean):
        threshold_range = [0.01, 0.1] + [i for i in range(1, 98)]
        qc_range = [self.wmean_qc_coef(wmean, threshold) for threshold in threshold_range]
        optimal_thr = [threshold_range[i] for i, j in enumerate(qc_range) if j == max(qc_range)][0]
        qc_coef = self.wmean_qc_coef(wmean, optimal_thr)
        return qc_coef

    def run_estimators(self):
        for name, estimator in self.dict_of_estimators.items():
            pu_clf = PUBaggingClassifier(name, estimator[0], estimator[1])
            pu_clf.run_bagging(self.X, self.y)
            prediction = pu_clf.get_df_with_qc(self.extended_list).fillna(-1)
            self.probas.append(prediction['Average_prob'].values)
            self.weights.append(prediction['QC_coef'].mean())
            self.thresholds.append(prediction['Optimal_threshold'].mean())
    
    def simple_weighted_mean(self, ind=True):
        weights = self.weights[ind]
        if np.sum(weights) == 0:
            weights = 1 * len(ind)
        return np.average(self.probas[ind], axis=0, weights=weights)
   
    def best_scored_proba(self):
        self.weights = np.array(self.weights)
        self.probas = np.array(self.probas)
        n_clfs = [i for i in range(len(self.probas))]
        max_score = 0
        comb_index = [list(combinations(n_clfs, i)) for i in range(1, len(n_clfs)+1)]
        comb_probs = [[self.probas[[*i]], [*i]] for comb_i in comb_index for i in comb_i]
        print('Finding best combination', end='\t')
        bar = Bar(max=len(comb_probs))
        print()
        for comb, ind in comb_probs:
            wmean = self.simple_weighted_mean(ind=ind)
            qc = self.perform_qc(wmean)
            if qc >= max_score:
                max_score = qc
                values = wmean
                index_max = ind
            bar.next()
        bar.finish()
        print(f'Best combination is: {np.array([key for key in MODELS.keys()])[index_max]}')
        print(f'QC: {max_score}')
        return pd.DataFrame({'gene_symbol':self.gene_names, 'Probability':values}).sort_values(by='Probability', ascending=False)
        
