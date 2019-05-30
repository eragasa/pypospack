import numpy as np
import pandas as pd
#https://towardsdatascience.com/machine-learning-with-python-easy-and-robust-method-to-fit-nonlinear-data-19e8a1ddbd49
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LassoCV
from sklearn.pipeline import make_pipeline

from pypospack.pyposmat.data import PyposmatDataFile

fn = "../../../data_test/Ni__eam__born_exp_bjs_00/data__Ni__eam__born_exp_bjs_04/pyposmat.kde.4.out"
data = PyposmatDataFile()
data.read(fn)

sqerr_names = ["{}.sqerr".format(q) for q in data.qoi_names]
df = data.df
df[sqerr_names] = np.square(data.df[data.error_names])

# Training set
test_set_fraction = 0.20
# Alpha (regularization strength) of LASSO regression
lasso_eps = 0.0001
lasso_alpha = 0.015
lasso_nalpha=20
lasso_tol=1e-6
lasso_iter=50000
# Min and max degree of polynomials features to consider
degree_min = 2
degree_max = 8
import numpy as np
nr,nc = df.shape

X_labels = list(sqerr_names)
X_labels.remove(X_labels[0])
y_label = data.error_names[0]
X=df[X_labels]
y=df[y_label]
# Test/train split
X_train, X_test, y_train, y_test = train_test_split(
        X,y,
        test_size=test_set_fraction)
# Make a pipeline model with polynomial transformation and LASSO regression with cross-validation, run it for increasing degree of polynomial (complexity of the model)

for degree in range(degree_min,degree_max+1):
    if degree == 2: 
        print('doing 2nd degree models')
    elif degree == 3: 
        print('doing 3rd degree models')
    else:
        print('doing {}th degree models.'.format(degree))
    model = make_pipeline(
            PolynomialFeatures(
                degree, 
                interaction_only=False), 
            LassoCV(
                eps=lasso_eps,
                fit_intercept=True,
                n_alphas=lasso_nalpha,
                max_iter=lasso_iter,
                normalize=True,cv=5)
            )
    model.fit(X_train,y_train)
    test_pred = np.array(model.predict(X_test))
    RMSE=np.sqrt(np.sum(np.square(test_pred-y_test)))
    test_score = model.score(X_test,y_test)
    print(test_pred)
    print(RMSE)
    print(test_score)

