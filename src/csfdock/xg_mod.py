from pandas import read_csv
from numpy import absolute
from matplotlib import pyplot
from numpy import mean
from numpy import std
from sklearn.datasets import make_regression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from xgboost import XGBRegressor
from sklearn.datasets import make_classification
from sklearn.model_selection import RepeatedStratifiedKFold
from xgboost import XGBClassifier
from matplotlib import pyplot


def ensamble(X, y):
    # define dataset
    X, y = make_regression(
        n_samples=1000, n_features=5, n_informative=5, noise=0.1, random_state=7
    )
    # define the model
    model = XGBRegressor()
    model.fit(X, y)
    # evaluate the model
    cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
    n_scores = cross_val_score(
        model,
        X,
        y,
        scoring="neg_mean_absolute_error",
        cv=cv,
        n_jobs=-1,
        error_score="raise",
    )
    # report performance
    print("MAE: %.3f (%.3f)" % (mean(n_scores), std(n_scores)))


def xg(file, ensemble=True, params=None):

    dataframe = read_csv(file, header=None)
    data = dataframe.values
    # split data into input and output columns
    X, y = data[1:, 1:-1], data[1:, -1]
    # define model
    model = XGBRegressor()
    model.fit(X, y)
    # define model evaluation method
    cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=10)
    # evaluate model
    scores = cross_val_score(
        model, X, y, scoring="neg_mean_absolute_error", cv=cv, n_jobs=-1
    )
    # force scores to be positive
    scores = absolute(scores)
    print("Mean MAE: %.3f (%.3f)" % (scores.mean(), scores.std()))

    if ensemble:
        ensamble(X, y)
    if params is not None:
        xg_param(X, y)


# explore xgboost number of trees effect on performance
def xg_param(X, y):
    def get_dataset():
        X, y = make_classification(
            n_samples=1000,
            n_features=5,
            n_informative=15,
            n_redundant=5,
            random_state=712,
        )
        return X, y

    # get a list of models to evaluate
    def get_models():
        trees = [10, 50, 100, 500, 1000, 5000]
        return {str(n): XGBClassifier(n_estimators=n) for n in trees}

    # evaluate a give model using cross-validation
    def evaluate_model(model):
        cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
        return cross_val_score(model, X, y, scoring="accuracy", cv=cv, n_jobs=-1)

    # define dataset
    X, y = get_dataset()
    # get the models to evaluate
    models = get_models()
    # evaluate the models and store results
    results, names = list(), list()
    for name, model in models.items():
        scores = evaluate_model(model)
        results.append(scores)
        names.append(name)
        print(">%s Accuracy: %.3f[mean] %.3f[std]" % (name, mean(scores), std(scores)))
    # plot model performance for comparison
    pyplot.boxplot(results, labels=names, showmeans=True)
    pyplot.show()
