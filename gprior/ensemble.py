from sklearn.calibration import CalibratedClassifierCV
from progress.bar import IncrementalBar, Bar
from itertools import combinations
import pandas as pd
import numpy as np
import random
import math
import warnings
import gprior.qctoolbox as QCtoolbox
import gprior.putoolbox as PUtoolbox
from gprior.var import *
warnings.filterwarnings("ignore", category=FutureWarning)



class PUBaggingClassifier():
    """
    PU-bagging classifier that implements and modify logic of Mordlet et al., 2014.
    We add hyperparameter tuning for each bootstrap iteration, implement custom cost
    function and expand list of classifiers.
    """
    def __init__(self, clf_name, clf_base, params, n_bootstrap, s_coef, tune, set_seed):
        self.clf_name = clf_name
        self.clf_base = clf_base
        self.params = params
        self.n_bootstrap = n_bootstrap
        self.set_seed = set_seed
        self.s_coef = s_coef
        self.tune = tune

    def run_bagging(self, X, y):
        """
        Kernel of PU-bagging classifier.
        Reference:
        Mordelet, F., Vert, J.-P. A bagging SVM to learn from positive and unlabeled examples. 
        Pattern Recognition Lett.(2013),http://dx.doi.org/10.1016/j.patrec.2013.06.010
        """
        self.gene_names = X.index 
        self.iP = y[y > 0].index
        self.iU = y[y <= 0].index
        num_oob = pd.DataFrame(np.zeros(shape=y.shape), index=y.index)
        sum_oob = pd.DataFrame(np.zeros(shape=y.shape), index=y.index)        
        print(f'{self.clf_name}', end='\t')
        bar = IncrementalBar(max=self.n_bootstrap)
        print()
        for num in range(self.n_bootstrap):
            if self.set_seed:
                np.random.seed(num)
            ib = np.random.choice(self.iU, replace=True, size=int(len(self.iP)*self.s_coef))
            i_oob = list(set(self.iU) - set(ib))        
            Xb = X[list(y > 0)].append(X.iloc[ib])
            yb = y[y > 0].append(y.iloc[ib])
            new_clf = self.clf_base
            
            if self.tune:
                new_clf = PUtoolbox.best_clf(self.clf_base,
                                        Xb,
                                        yb,
                                        PARAMS[self.clf_name])
            
            if self.set_seed:
                new_clf.set_params(random_state=num)
                
            new_clf = CalibratedClassifierCV(new_clf, cv=3)
            new_clf.fit(Xb, yb)
            sum_oob.loc[i_oob, 0] += new_clf.predict_proba(X.iloc[i_oob])[:,1]
            num_oob.loc[i_oob, 0] += 1
            bar.next()

        probs = (100 * sum_oob / num_oob).fillna(-1)[0]
        bar.finish()
        return pd.DataFrame({'gene_symbol': self.gene_names,
                               'average_prob': probs})
        
class EnsembleClassifier():    
    """
    Ensemble class that aggrigate pu-bagging classification results from 
    different classification algorithms.
    Perform qc and 'optimal combination' finding strategy.
    """
    def __init__(self, X, y, dict_of_estimators, set_seed=False, tune=False, n_bootstrap=10, s_coef=1):
        self.X = X
        self.y = y
        self.set_seed = set_seed
        self.dict_of_estimators = dict_of_estimators
        self.n_bootstrap = n_bootstrap
        self.alg_eval_set = False
        self.s_coef = s_coef
        self.tune = tune
    
    def set_aes(self, alg_eval_set):
        self.alg_eval_set = set(alg_eval_set.gene_symbol.values) 

    def set_ytrue(self):
        self.true_y = QCtoolbox.give_true_y(self.X, self.y, self.alg_eval_set)

    def run_estimators(self):
        """
        Iteratively run PU-bagging with all provided ML algorithms,
        save predictions and weights (PU-scores).
        """
        
        self.weights = list()
        self.probas = list()
        for name, estimator in self.dict_of_estimators.items():
            pu_clf = PUBaggingClassifier(name,
                                         estimator[0],
                                         estimator[1],
                                         set_seed=self.set_seed,
                                         tune=self.tune,
                                         n_bootstrap=self.n_bootstrap,
                                         s_coef=self.s_coef)
            
            prediction = pu_clf.run_bagging(self.X, self.y).average_prob.values
            
            if self.alg_eval_set:
                weight = QCtoolbox.give_summary(wmean=prediction,
                                            threshold_range=THRESHOLD_RANGE,
                                            true_y=self.true_y)
                self.weights.append(weight)
                print(f'PU-score: {weight}', end='\n\n')
            
            self.probas.append(prediction)        
    
    def prediction_without_aes(self):
        self.probas = np.array(self.probas)
        n_clfs = list(range(len(self.dict_of_estimators)))
        return pd.DataFrame({
                'gene_symbol':self.X.index,
                'LR': self.probas[0],
                'SVM': self.probas[1],
                'Ada': self.probas[2],
                'RF': self.probas[3],
                'DT': self.probas[4],
                'wmean': PUtoolbox.simple_weighted_mean(ind=n_clfs,
                                                        probas=self.probas)
            }).sort_values(by='wmean', ascending=False)

    def prediction_with_aes(self):
        """
        Find optimal combination of weighted predictions.
        """
        self.weights = np.array(self.weights)
        self.probas = np.array(self.probas)
        n_clfs = list(range(len(self.dict_of_estimators)))
        power_set = PUtoolbox.powerset(n_clfs)
        max_score, values, index_max = PUtoolbox.give_max_score(power_set, self.true_y, self.probas)
        print(f'Best combination is: {np.array(list(MODELS.keys()))[index_max]}')
        print(f'PU-score: {max_score}')
        return pd.DataFrame({'gene_symbol':self.X.index,
                             'Probability':values}) \
                             .sort_values(by='Probability',
                                          ascending=False)