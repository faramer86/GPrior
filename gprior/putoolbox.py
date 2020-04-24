from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer
from itertools import combinations
from gprior.qctoolbox import cost_function

def best_clf(model, X_train, y_train, param_dist, cv=3):
    """
    Hyperparameter tuning based on grid search.
    It use PU_score as scoring function.
    """
    clf = GridSearchCV(model,
                       param_dist,
                       cv=cv,
                       scoring = make_scorer(cost_function,
                                             greater_is_better=True),
                       iid=False,
                       n_jobs=-1)
    clf.fit(X_train, y_train)
    return clf.best_estimator_

def powerset(someset):
    size = len(someset)
    return list(list(sub) for k in range(1, size+1) \
            for sub in combinations(someset,k))[1:]