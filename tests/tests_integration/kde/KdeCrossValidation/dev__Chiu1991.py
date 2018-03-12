# Reference: Chiu 1991
# https://www.r-bloggers.com/cross-validation-for-kernel-density-estimation/
import numpy as np
from scipy import stats
from scipy import integrate
from scipy import optimize

#def Silverman1986_h(X):
#    std = np.std(X)
#    X_shape = X.shape
#    N=X_shape[0]
#    return 1.06*std*N**(-.2)

def Silverman1986_h(X):
    kde = stats.gaussian_kde(X,'silverman')
    return kde.factor

def Scott1986_h(X):
    kde = stats.gaussian_kde(X,'scott')

def Chiu1999_crossvalidation(X):
     
    def fhati(X,h,i):
        if h is float: _h=h
        else: _h=h[0]

        Xi = np.delete(X,i)
        kde = stats.gaussian_kde(Xi,_h)
        return kde(X[i])

    def J(X,h):
        if h is float: _h=h
        else: _h=h[0]

        fhat = stats.gaussian_kde(X,_h)
        #F1 = integrate.quad(lambda x: fhat(x)**2,-np.inf,np.inf)[0]
        F1 = fhat.integrate_kde(fhat)
        F2 = np.array([fhati(X,h,i) for i in range(X.shape[0])])
        return F1-2*np.mean(F2)

    #h0 = Silverman1986_h(X)
    h0 = .5
    results = optimize.minimize(
        lambda h: J(X,h),
        h0,
        method='Nelder-Mead')

    return results.x[0]
    
if __name__ == "__main__":
    #mu1, sigma1, N1 = 0, 1, 50
    #mu2, sigma2, N2 = 1, 1, 50

    #X1 = np.random.normal(mu1,sigma1,50)
    #X2 = np.random.normal(mu2,sigma2,50)
    #X = np.concatenate((X1,X2))

    from pypospack.pyposmat.data import PyposmatDataFile

    for i in range(9):
        f = '../../../data_test/Ni__eam__born_exp_fs_00/data__Ni__eam__born_exp_fs_02/pyposmat.kde.{}.out'.format(i) 
        data=PyposmatDataFile()
        data.read(filename=fn)
        X=data.df[data.parameter_names].values

        print("Silverman1986_h:{}".format(Silverman1986_h(X.T)))
        print("Scott1986_h:{}".format(Scott1986_h(X.T)))
        h = Chiu1999_crossvalidation(X.T)
        print('Chiu1999_h:{}'.format(h))
        print(80*'-')
