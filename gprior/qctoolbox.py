import numpy as np
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from gprior.var import *

def pu_score(true_y, pred_y):
    """
    Quality score that is used for score function.
    Reference:
    Wee Sun Lee and Bing Liu. Learning with positive and unlabeled exam-
    ples using weighted logistic regression. In Proceedings of the Twentieth
    International Conference on Machine Learning (ICML), pages 448–455, 2003.
    """
    Pr = np.sum(pred_y) / pred_y.size
    recall = recall_score(true_y, pred_y)
    return recall**2/Pr

def cost_function(true_y, pred_y):
    """
    Custom cost function for hyperparametrs tuning.
    It based on PU_score = recall^2/Pr[y==1]
    See pu_score function for references.
    """
    true_y, pred_y = np.array(true_y), np.array(pred_y)
    if np.sum(pred_y) == 0:
        return 0
    return pu_score(true_y, pred_y)

def give_true_y(X, y, alg_eval_set):
    """
    Return np.array with true answers for the whole geneset
    based on algorithm selection set.
    Used for Pr[y==1] estimation.
    See original article for more details.
    """
    gene_names = X.index 
    iP = y[y > 0].index
    true_y = list(map(lambda gene: (gene in alg_eval_set) & \
        (gene not in gene_names[iP]), gene_names))
    return np.array(true_y) 

def wmean_qc_coef(wmean, threshold, true_y):
    """
    Return PU score for different threshold values.
    Used for optimal threshold estimation.
    """
    pred_y = np.array(wmean >= threshold)
    if np.sum(pred_y) == 0:
        return 0
    return pu_score(true_y, pred_y)
    
def give_qc_range(wmean, threshold_range, true_y):
    """
    Return list with pu-scores for different threshold values.
    """
    return [wmean_qc_coef(wmean, threshold, true_y)
            for threshold in threshold_range]

def give_summary(wmean, true_y, thr=False, eAUC=False):
    """
    Depends on thr value, it return maximal PU score
    or maximal PU score with optimal probability threshold
    """
    threshold_range = wmean[::3]
    qc_range = give_qc_range(wmean, threshold_range, true_y)
    qc_coef = max(qc_range)
    if eAUC:
        qc_coef = roc_auc_score(true_y, wmean)
        return qc_coef
    if thr:
        optimal_thr = threshold_range[qc_range.index(qc_coef)] 
        return optimal_thr, qc_coef
    return qc_coef

