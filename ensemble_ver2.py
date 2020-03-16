from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer
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
import qctoolbox
warnings.filterwarnings("ignore", category=FutureWarning)

class PUBaggingClassifier():
    """
    PU-bagging classifier that implements and modify logic of Mordlet et al., 2014.
    We add hyperparameter tuning for each bootstrap iteration, implement custom cost
    function and expand list of classifiers.
    """
    def __init__(self, clf_name, clf_base, params, n_bootstrap, all_genes):
        self.clf_name = clf_name
        self.clf_base = clf_base
        self.params = params
        self.n_bootstrap = n_bootstrap
        self.all_genes = all_genes

    def best_clf(self, model, X_train, y_train, param_dist):
        """
        Hyperparameter tuning based on grid search.
        It use PU_score as scoring function.
        """
        clf = GridSearchCV(model,
                           param_dist,
                           cv=2,
                           scoring = make_scorer(qctoolbox.cost_function,
                                                  greater_is_better=True),
                           iid=False,
                           n_jobs=-1)
        clf.fit(X_train, y_train)
        return clf.best_estimator_

    def run_bagging(self, X, y):
        """
        Kernel of PU-bagging classifier.
        Reference:
        Mordelet, F., Vert, J.-P. A bagging SVM to learn from positive and unlabeled examples. 
        Pattern Recognition Lett.(2013),http://dx.doi.org/10.1016/j.patrec.2013.06.010
        """
        y = pd.Series(list(gtex_frames['gene_symbol'].map(lambda x: 1 if x in scz_genes else 0)))
        self.gene_names = X.index 
        self.iP = y[y > 0].index
        self.iU = y[y <= 0].index        
        num_oob = pd.DataFrame(np.zeros(shape=y.shape), index=y.index)
        sum_oob = pd.DataFrame(np.zeros(shape=y.shape), index=y.index)        
        print(f'{self.clf_name}', end='\t')
        bar = IncrementalBar(max=self.n_bootstrap)
        print()
        for num in range(self.n_bootstrap):
            np.random.seed(num)
            
            ib = np.random.choice(self.iU, replace=True, size=len(self.iP))
            i_oob = list(set(self.iU) - set(ib))        
            
            # pib = np.random.choice(self.iP, replace=False, size=int(len(self.iP)*0.5))
            # pi_oob = list(set(self.iP) - set(pib))
            
            Xb = X[list(y>0)].append(X.iloc[ib])
            yb = y[y > 0].append(y.iloc[ib])
            new_clf = self.clf_base
            # new_clf = self.best_clf(self.clf_base,
            #                         Xb,
            #                         yb,
            #                         PARAMS[self.clf_name])
            
            new_clf.set_params(random_state=num)            
            
            new_clf.fit(Xb, yb)
            #sum_oob.loc[i_oob, 0] += new_clf.predict_proba(X.iloc[i_oob])[:,1]
            #num_oob.loc[i_oob, 0] += 1
            sum_oob.loc[i_oob, num] = new_clf.predict_proba(X.iloc[i_oob])[:,1]
            bar.next()
            #print(f'{num}', end=' ')

        df = sum_oob.set_index(self.gene_names)
        weights = [qctoolbox.give_summary(df[i], THRESHOLD_RANGE, self.algeval) for i in range(len(list(df)))]
        print(weights)
        ind = [ind for ind, w in enumerate(weights) if w > 1]
        wmean = (100 * df[ind] * np.array(weights)[ind]).sum(axis=1)/sum(weights)
        mean = 100*df[ind].sum(axis=1)/len(ind)
        print(mean)
        #probs = (100 * sum_oob / num_oob).fillna(-1)[0]
        clf_df = pd.DataFrame({'gene_symbol': self.gene_names,
                               'average_prob': list(mean)}) # probs
        bar.finish()
        return clf_df
        
class EnsembleClassifier():    
    """
    Ensemble class that aggrigate pu-bagging classification results from 
    different classification algorithms.
    Perform qc and 'optimal combination' finding strategy.
    """
    def __init__(self, X, y, dict_of_estimators, alg_eval_set, true_genes, n_bootstrap=100):
        self.X = X
        self.y = y
        self.alg_eval_set = set(alg_eval_set.gene_symbol.values)
        self.comb = np.concatenate([list(alg_eval_set), list(true_genes)])
        self.true_y = qctoolbox.give_true_y(self.X, self.y, self.alg_eval_set)
        self.dict_of_estimators = dict_of_estimators
        self.n_bootstrap = n_bootstrap
        
    def run_estimators(self):
        """
        Iteratively run PU-bagging with all provided ML algorithms,
        save predictions and weights (PU-scores).
        """
        self.probas = list()
        self.weights = list()
        for name, estimator in self.dict_of_estimators.items():
            pu_clf = PUBaggingClassifier(name,
                                         estimator[0],
                                         estimator[1],
                                         n_bootstrap=self.n_bootstrap,
                                         all_genes = self.comb) #algeval = self.true_y
            prediction = pu_clf.run_bagging(self.X, self.y).average_prob.values
            weight = qctoolbox.give_summary(prediction, THRESHOLD_RANGE, self.true_y)
            self.probas.append(prediction)
            self.weights.append(weight)
            print(f'PU-score: {weight}', end='\n\n')
    
    def simple_weighted_mean(self, ind):
        """
        Calculate weighted mean for each combination.
        """
        # weights = self.weights[ind]
        # if np.sum(weights) == 0:
        #     weights = 1 * len(ind)
        return np.average(self.probas[ind], axis=0)  #, weights=weights
   
    def best_scored_proba(self):
        """
        Find optimal combination of weighted predictions.
        """
        self.weights = np.array(self.weights)
        self.probas = np.array(self.probas)
        n_clfs = [i for i in range(len(self.probas))]
        max_score = 0
        comb_index = [list(combinations(n_clfs, i)) for i in range(1, len(n_clfs)+1)]
        comb_probs = [[self.probas[[*i]], [*i]] for comb_i in comb_index for i in comb_i]
        print('Finding best combination', end='\t')
        bar = Bar(max=len(comb_probs))
        print()
        for comb, ind in comb_probs: # comb?
            wmean = self.simple_weighted_mean(ind=ind)
            thr, qc = qctoolbox.give_summary(wmean, THRESHOLD_RANGE, self.true_y, thr=True)
            if qc >= max_score:
                max_score = qc
                values = wmean
                index_max = ind
                best_thr = thr
            bar.next()
        bar.finish()
        print(f'Best combination is: {np.array(list(MODELS.keys()))[index_max]}')
        print(f'PU-score: {max_score}')
        #print(f'Optimal threshold: {best_thr}')
        return pd.DataFrame({'gene_symbol':self.X.index,
                             'Probability':values}).sort_values(by='Probability', ascending=False)


if __name__ == "__main__":
        
    import numpy as np
    from itertools import combinations
    from sklearn.cluster import FeatureAgglomeration
    from sklearn.preprocessing import RobustScaler
    import pandas as pd

    def return_x_y(gtex_frames, scz_genes):
        y = pd.Series(list(gtex_frames['gene_symbol'].map(lambda x: 1 if x in scz_genes else 0)))
        #gtex_frames.drop('gene_id', axis=1, inplace=True)
        x = gtex_frames.fillna(0).set_index('gene_symbol')
        return x, y

    gtex_frames = pd.read_csv('/home/nikita/Desktop/triple_phe/SCZ/new_combined_scz.tsv',sep='\t')
    true_genes = pd.read_csv('/home/nikita/Documents/work/git_projects/GWASPriors/SCZ/SCZ/true_genes.csv')
    scz_genes = set(true_genes[true_genes['Gene(s)'].map(lambda x: ',' not in str(x))]['Gene(s)']) & set(gtex_frames['gene_symbol'])
    gtex_frames = gtex_frames.loc[gtex_frames.apply(lambda x: sum(x == 0), axis=1) <= 80,]
    extended = pd.read_csv('/home/nikita/Documents/work/git_projects/GWASPriors/SCZ/SCZ/extended.csv')
    ext_genes = extended['GWAS_mapper_predictions_PPI']
    #ext_genes = pd.read_csv('/home/nikita/Desktop/illustration/ppi.tsv')
    gsea = set(ext_genes.values)
    gsea = pd.DataFrame({'gene_symbol':pd.Series(list(gsea))})
    x, y = return_x_y(gtex_frames, scz_genes)
    gene_names = x.index.copy()       
    scaler = RobustScaler()
    x = pd.DataFrame(FeatureAgglomeration(n_clusters=25).fit_transform(scaler.fit_transform(x)), index=gene_names)
    ens = EnsembleClassifier(x, y, MODELS, gsea, scz_genes, n_bootstrap=60)
    ens.run_estimators()
    ens.best_scored_proba()