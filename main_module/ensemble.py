from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import make_scorer
from sklearn.model_selection import RandomizedSearchCV
from sklearn.datasets import make_blobs
from progress.bar import IncrementalBar, Bar
from sklearn.base import clone
from itertools import combinations
import pandas as pd
import numpy as np
import random
import math
import warnings
import main_module.qctoolbox as qctoolbox
from main_module.var import *
warnings.filterwarnings("ignore", category=FutureWarning)

class PUBaggingClassifier():
    """
    PU-bagging classifier that implements and modify logic of Mordlet et al., 2014.
    We add hyperparameter tuning for each bootstrap iteration, implement custom cost
    function and expand list of classifiers.
    """
    def __init__(self, clf_name, clf_base, params, n_bootstrap, set_seed):
        self.clf_name = clf_name
        self.clf_base = clf_base
        self.params = params
        self.n_bootstrap = n_bootstrap
        self.set_seed = set_seed

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
            ib = np.random.choice(self.iU, replace=True, size=len(self.iP))
            i_oob = list(set(self.iU) - set(ib))        
            Xb = X[list(y > 0)].append(X.iloc[ib])
            yb = y[y > 0].append(y.iloc[ib])
            new_clf = self.clf_base
            
            # new_clf = self.best_clf(self.clf_base,
            #                         Xb,
            #                         yb,
            #                         PARAMS[self.clf_name])
            
            if self.set_seed:
                new_clf.set_params(random_state=num)
            new_clf = CalibratedClassifierCV(new_clf, cv=3)
            new_clf.fit(Xb, yb)
            sum_oob.loc[i_oob, 0] += new_clf.predict_proba(X.iloc[i_oob])[:,1]
            num_oob.loc[i_oob, 0] += 1
            bar.next()
        probs = (100 * sum_oob / num_oob).fillna(-1)[0]
        clf_df = pd.DataFrame({'gene_symbol': self.gene_names,
                               'average_prob': probs})
        bar.finish()
        return clf_df
        
class EnsembleClassifier():    
    """
    Ensemble class that aggrigate pu-bagging classification results from 
    different classification algorithms.
    Perform qc and 'optimal combination' finding strategy.
    """
    def __init__(self, X, y, dict_of_estimators, alg_eval_set, set_seed=False, n_bootstrap=10):
        self.X = X
        self.y = y
        self.set_seed = set_seed
        self.alg_eval_set = set(alg_eval_set.gene_symbol.values)
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
                                         set_seed=self.set_seed,
                                         n_bootstrap=self.n_bootstrap)
            prediction = pu_clf.run_bagging(self.X, self.y).average_prob.values
            weight = qctoolbox.give_summary(wmean=prediction,
                                            threshold_range=THRESHOLD_RANGE,
                                            true_y=self.true_y)
            self.probas.append(prediction)
            self.weights.append(weight)
            print(f'PU-score: {weight}', end='\n\n')
    
    def simple_weighted_mean(self, ind, use_weights=False):
        """
        Calculate weighted mean for each combination.
        """
        if use_weights:
            weights = self.weights[ind]
            if np.sum(weights) == 0:
                weights = [1] * len(ind)
            return np.average(self.probas[ind],
                              axis=0,
                              weights=weights)  #, weights=weights
        return np.average(self.probas[ind], axis=0)

    def best_scored_proba(self):
        """
        Find optimal combination of weighted predictions.
        """
        self.weights = np.array(self.weights)
        self.probas = np.array(self.probas)
        n_clfs = [i for i in range(len(self.probas))]
        index_sets = [list(combinations(n_clfs, i)) for i in range(1, len(n_clfs)+1)]
        comb_probs = [[*i] for comb in index_sets for i in comb]
        print('Finding best combination', end='\t')
        bar = Bar(max=len(comb_probs))
        print()
        max_score = 0
        for ind in comb_probs:
            wmean = self.simple_weighted_mean(ind=ind)
            qc = qctoolbox.give_summary(wmean=wmean,
                                        threshold_range=THRESHOLD_RANGE,
                                        true_y=self.true_y,
                                        thr=False)
            if qc >= max_score:
                max_score = qc
                values = wmean
                index_max = ind
                #best_thr = thr
            bar.next()
        bar.finish()
        print(f'Best combination is: {np.array(list(MODELS.keys()))[index_max]}')
        print(f'PU-score: {max_score}')
        #print(f'Optimal threshold: {best_thr}')
        return pd.DataFrame({'gene_symbol':self.X.index,
                             'Probability':values}) \
                             .sort_values(by='Probability',
                                          ascending=False)


if __name__ == "__main__":
    pass
        
    # import numpy as np
    # from itertools import combinations
    # from sklearn.cluster import FeatureAgglomeration
    # from sklearn.preprocessing import RobustScaler
    # import pandas as pd

    # def return_x_y(gtex_frames, scz_genes):
    #     y = pd.Series(list(gtex_frames['gene_symbol'].map(lambda x: 1 if x in scz_genes else 0)))
    #     #gtex_frames.drop('gene_id', axis=1, inplace=True)
    #     x = gtex_frames.fillna(0).set_index('gene_symbol')
    #     return x, y

    # gtex_frames = pd.read_csv('/home/nikita/Desktop/triple_phe/SCZ/new_combined_scz.tsv',sep='\t')
    # true_genes = pd.read_csv('/home/nikita/Documents/work/git_projects/GWASPriors/SCZ/SCZ/true_genes.csv')
    # scz_genes = set(true_genes[true_genes['Gene(s)'].map(lambda x: ',' not in str(x))]['Gene(s)']) & set(gtex_frames['gene_symbol'])
    # gtex_frames = gtex_frames.loc[gtex_frames.apply(lambda x: sum(x == 0), axis=1) <= 80,]
    # extended = pd.read_csv('/home/nikita/Documents/work/git_projects/GWASPriors/SCZ/SCZ/extended.csv')
    # ext_genes = extended['GWAS_mapper_predictions_PPI']
    # #ext_genes = pd.read_csv('/home/nikita/Desktop/illustration/ppi.tsv')
    # gsea = set(ext_genes.values)
    # gsea = pd.DataFrame({'gene_symbol':pd.Series(list(gsea))})
    # x, y = return_x_y(gtex_frames, scz_genes)
    # gene_names = x.index.copy()       
    # scaler = RobustScaler()
    # x = pd.DataFrame(FeatureAgglomeration(n_clusters=25).fit_transform(scaler.fit_transform(x)), index=gene_names)
    # ens = EnsembleClassifier(x, y, MODELS, gsea, n_bootstrap=60)
    # ens.run_estimators()
    # ens.best_scored_proba()
    def wmean_qc_coef(wmean, threshold, extended_list, true_genes):
        wmean_names = wmean['gene_symbol'].values
        pred_y = np.array(wmean['Probability'].values >= threshold)
        new_y = np.array(list(map(lambda gene: (gene in extended_list['gene_symbol'].values) & \
                            (gene not in true_genes['gene_symbol'].values), wmean_names))) #wmean_names
        
        if np.sum(pred_y) == 0:
            return 0
        Pr = np.sum(pred_y) / pred_y.size
        TP = np.sum(new_y & pred_y) # TP
        FN = np.sum(new_y) - TP # FN
        recall = TP/(TP+FN)
        coef = recall**2/Pr
        norm_coef = coef * np.sum(new_y) / new_y.size
        return coef # coef
    
    def perform_qc(wmean, extended_list, true_genes):
        threshold_range = wmean.Probability.values #[0.01, 0.1] + [i for i in range(1, 98)]
        qc_range = [wmean_qc_coef(wmean, threshold, extended_list, true_genes) for threshold in threshold_range]
        #print(qc_range)
        optimal_thr = [threshold_range[i] for i, j in enumerate(qc_range) if j == max(qc_range)][0]
        qc_coef = wmean_qc_coef(wmean, optimal_thr, extended_list, true_genes)
        print('qc_coef:', qc_coef)
        print('opt_threshols:', optimal_thr)
        return qc_range

    def plot_precision_recall_vs_threshold(precisions, recalls, thresholds):
        """
        Modified from:
        Hands-On Machine learning with Scikit-Learn
        and TensorFlow; p.89
        """
        plt.figure(figsize=(3, 3))
        plt.title("Precision and Recall Scores as a function of the decision threshold")
        plt.plot(thresholds, precisions[:-1], "b--", label="Precision")
        plt.plot(thresholds, recalls[:-1], "g-", label="Recall")
        plt.ylabel("Score")
        plt.xlabel("Decision Threshold")
        plt.legend(loc='best')
        plt.savefig('ProbThresh_PU.pdf')

    from sklearn.cluster import FeatureAgglomeration
    from sklearn.preprocessing import RobustScaler
    import pandas as pd
    from ensemble import *

    ntrain = [round(212*dob/100) for dob in range(2,22,2)]

    import math
    from sklearn.preprocessing import StandardScaler

    roc_dict = dict()
    f_dict = dict()
    pu_dict = dict()

    X = pd.read_csv('/home/nikita/Downloads/cancer.csv')
    y = X.diagnosis
    del X['diagnosis']
    del X['Unnamed: 32']
    del X['id']
    scaler = StandardScaler()
    scaler.fit(X)
    X = pd.DataFrame(scaler.transform(X))
    y = y.map(lambda x: 1 if x == 'M' else 0)

    for ntr in  ntrain:
        
        roc_store = list()
        f1_store = list()
        pu_store = list()
        hidden_size = y.sum() - ntr

        print(f'NTR ==== {ntr}')
        for i in range(100):
            
            y_orig = y.copy()
            y_copy = y.copy()

            y_copy.loc[
                    np.random.choice(
                        y_copy[y_copy == 1].index, 
                        replace = False, 
                        size = hidden_size
                    )
                ] = 0
            val_y = y_copy[y_copy == 1].copy()
            y_copy.loc[
                    np.random.choice(
                        y_copy[y_copy == 1].index, 
                        replace = False, 
                        size = round(ntr*1/2)
                    )
                ] = 0
            als = pd.DataFrame({'gene_symbol':list(set(val_y[val_y==1].index)-set(y_copy[y_copy==1].index))})
            X['gene_symbol'] = X.index
            ens = EnsembleClassifier(X, y_copy, MODELS, als, n_bootstrap=25)
            ens.run_estimators()
            wmean = ens.best_scored_proba()
            wmean['Probability'] = wmean.Probability.apply(lambda x: 100 if x==-1 else x)
            extended_list = pd.DataFrame({'gene_symbol':y_orig[y_orig==1].index.values})
            true_genes = pd.DataFrame({'gene_symbol':val_y.index.values})
            qc = max(perform_qc(wmean, extended_list, true_genes))
            from sklearn.metrics import *

            from math import isclose
            st_res = pd.DataFrame({'prob':wmean.Probability}).dropna()
            p, r, thresholds = precision_recall_curve(y_orig.values[st_res.index], st_res)
            opt_threshold = thresholds[list(map(lambda n: isclose(n[0], n[1], rel_tol=0.01), zip(p,r)))[:-1]]
            roc_auc = roc_auc_score(y_orig.values[st_res.index], st_res>=opt_threshold[0]) 
            f1 = f1_score(y_orig.values[st_res.index], st_res>=opt_threshold[0])
            pu_score = recall_score(y_orig.values[st_res.index], st_res>=opt_threshold[0])**2/(int((st_res>=opt_threshold[0]).sum())/len(st_res))
            print(f'ALT PU = {pu_score}' )

            pu_store.append(qc)
        pu_dict[str(f'{ntr}')] = pu_store