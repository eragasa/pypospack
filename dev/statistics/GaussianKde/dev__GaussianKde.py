import numpy as np
from scipy import stats
from pypospack.statistics import GaussianKde

def measure(n):
    "Measurement model, return two coupled measurements."
    m1 = np.random.normal(size=n)
    m2 = np.random.normal(scale=0.5, size=n)
    X = np.vstack([m1+m2, m1-m2])
    return X
if __name__ == "__main__":

    X = measure(2000)
    print('type(X):',type(X).__name__)
    print('X.shape:',X.shape)
    scipy_kde = stats.gaussian_kde(X)
    print('scipy_kde.d:',scipy_kde.d)
    print('scipy_kde.n:',scipy_kde.n)
    pypospack_kde = GaussianKde(X)
    print('pypospack_kde.d:',pypospack_kde.d)
    print('pypospack_kde.n:',pypospack_kde.n)


    print('----test mean of standard scaler----')
    print(pypospack_kde.scaler.mean_)
    print(X[0,:].mean()-pypospack_kde.scaler.mean_[0])
    print(X[1,:].mean()-pypospack_kde.scaler.mean_[1])

    print('----- test variance of standard scaler----')
    print(np.sqrt(pypospack_kde.scaler.var_))
    print(X[0,:].std()-np.sqrt(pypospack_kde.scaler.var_[0]))
    print(X[1,:].std()-np.sqrt(pypospack_kde.scaler.var_[1]))
   

    print('----- test evaluation of points -----')
    points = measure(10)
    print(points.shape)

    P_of_X_scipy = scipy_kde.evaluate(points)
    P_of_X_pypospack = pypospack_kde.evaluate(points)
    
    print(P_of_X_scipy)
    print(P_of_X_pypospack)
    diff = scipy_kde.evaluate(points)-pypospack_kde.evaluate(points)
    print([k for k in diff.tolist()])

    print('----- resample points ---------------')
    n = 5
    sci_points = scipy_kde.resample(n)
    pyp_points = pypospack_kde.resample(n)

    print(sci_points.shape)
    print(pyp_points.shape)
    assert sci_points.shape == pyp_points.shape

    print(sci_points[0,:].mean() - sci_points[0,:].mean())
    print(sci_points[1,:].mean() - pyp_points[1,:].mean())
