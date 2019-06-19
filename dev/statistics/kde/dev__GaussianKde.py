import numpy as np
from scipy import stats
from pypospack.statistics import GaussianKde

def measure(n):
    "Measurement model, return two coupled measurements."
    m1 = np.random.normal(size=n)
    m2 = np.random.normal(scale=0.5, size=n)
    return m1+m2, m1-m2

if __name__ == "__main__":

    m1,m2 = measure(2000)
    values = np.vstack([m1,m2])
    
    X = values 

    scipy_kde = stats.gaussian_kde(values)
    pypospack_kde = GaussianKde(values)

    print('----test mean of standard scaler----')
    print([m1.mean(),m2.mean()])
    print(pypospack_kde.scaler.mean_)
    print(m1.mean()-pypospack_kde.scaler.mean_[0])
    print(m2.mean()-pypospack_kde.scaler.mean_[1])

    print('----- test variance of standard scaler----')
    print([m1.std(),m2.std()])
    print(np.sqrt(pypospack_kde.scaler.var_))
    print(m1.std()-np.sqrt(pypospack_kde.scaler.var_[0]))
    print(m2.std()-np.sqrt(pypospack_kde.scaler.var_[1]))
   

    print('----- test evaluation of points -----')
    points = np.array([2,2])
    print(scipy_kde.evaluate(points)-pypospack_kde.evaluate(points))

    print('----- resample points ---------------')
    sci_points = scipy_kde.resample(20000)
    pyp_points = pypospack_kde.resample(20000)

    print(sci_points.shape)
    print(pyp_points.shape)

    print(sci_points[0,:].mean() - sci_points[0,:].mean())
    print(sci_points[1,:].mean() - pyp_points[1,:].mean())
