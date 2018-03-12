import numpy as np
from scipy import stats
from scipy import integrate
from scipy import optimize

def Silverman1986_h(X):
    """
    Get smoothing factor using Silverman method
    """
    kde = stats.gaussian_kde(X,'silverman')
    return kde.factor

def Chiu1999_h(X):
    """
    Cross validation method of Chiu 1999

    Chiu, S.T.,  Annls. of Stat., 1991, 19, 1883-1905
    https://projecteuclid.org/download/pdf_1/euclid.aos/1176348376
    """
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
