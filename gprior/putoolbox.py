from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer
from itertools import combinations
import gprior.qctoolbox as QCtoolbox
from progress.bar import Bar
import numpy as np

def best_clf(model, X_train, y_train, param_dist, cv=3):
    """
    Hyperparameter tuning based on grid search.
    It use PU_score as scoring function.
    """
    clf = GridSearchCV(model,
                       param_dist,
                       cv=cv,
                       scoring = make_scorer(QCtoolbox.cost_function,
                                             greater_is_better=True),
                       iid=False,
                       n_jobs=-1)
    clf.fit(X_train, y_train)
    return clf.best_estimator_

def powerset(someset):
    size = len(someset)
    return list(list(sub) for k in range(1, size+1) \
            for sub in combinations(someset,k))

def simple_weighted_mean(ind, probas, weights=None, use_weights=False):
    """
    Calculate weighted mean for each combination.
    """
    if use_weights:
        weights = weights[ind]
        if np.sum(weights) == 0:
            weights = [1] * len(ind)
        return np.average(probas[ind],
                          axis=0,
                          weights=weights)
    return np.average(probas[ind], axis=0)

def give_max_score(power_set, true_y, probas, weights=None):

    print('Finding best combination', end='\t')
    bar = Bar(max=len(power_set))
    print()
        
    max_score = 0
    for ind in power_set:
        wmean = simple_weighted_mean(ind, probas, weights)
        qc = QCtoolbox.give_summary(wmean=wmean,
                                    true_y=true_y,
                                    thr=False)
        if qc >= max_score:
            max_score = qc
            values = wmean
            index_max = ind
        bar.next()
    bar.finish()
    return max_score, values, index_max
