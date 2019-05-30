"""
example on how to GMM

adapted from
Jake VanderPlass, Python Data Science Handbook
https://jakevdp.github.io/PythonDataScienceHandbook/05.12-gaussian-mixtures.html
"""
import os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns; sns.set()
import numpy as np
import pandas as pd


from sklearn.mixture import GMM
from sklearn.mixture import GaussianMixture
def plot_ellipse(position,covariance,ax=None,**kwargs):
    if ax is None:
        fig, ax = plt.subplots(1,1)

    if covariance.shape == (2,2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1,0],U[0,0]))
        width, height = 2 * np.sqrt(s)
    else:
        angle = 0
        width, height = 2 * np.sqrt(covariance)

    # draw ellipse

    for nsig in range(1,4):
        ax.add_patch(Ellipse(position,nsig*width,nsig*height,angle,**kwargs))

def plot_gmm(gmm,X,label=True,ax=None):
    if ax is None:
        fig, ax = plt.subplots(1,1)

    cluster_id = gmm.fit(X).predict(X)

    x = X[:,0]
    y = X[:,1]
    if label:
        ax.scatter(x,y,c=cluster_id,s=40,cmap='viridis',zorder=2)
    else:
        ax.scatter(x,y,s=40,zorder=2)

    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_,gmm.covariances_,gmm.weights_):
        plot_ellipse(pos,covar,alpha=w*w_factor,ax=ax)
    plt.show()

if __name__ == "__main__":
    # generate data
    from sklearn.datasets.samples_generator import make_blobs
    X, y_true = make_blobs(n_samples=400,centers=4,cluster_std=0.60,random_state=0)

    print('X.shape:{}'.format(X.shape))
    X = X[:,::-1]
    print('X.shape:{}'.format(X.shape))
    gmm = GaussianMixture(n_components=4).fit(X)
    cluster_id = gmm.predict(X)

    fig, ax = plt.subplots(1,1)
    x = X[:,0]
    y = X[:,1]
    ax.scatter(x,y,c=cluster_id,s=40,cmap='viridis')
    plt.show()
    plt.close()

    probs = gmm.predict_proba(X)
    print(probs[:5].round(3))
    n_samples, n_clusters = probs.shape
    print('n_samples={}'.format(n_samples))
    print('n_clusters={}'.format(n_clusters))

    gmm = GaussianMixture(n_components=4, covariance_type='full', random_state=42)
    rng = np.random.RandomState(13)
    X_stretched = np.dot(X, rng.randn(2, 2))
    plot_gmm(gmm, X_stretched)

    print('make_moons test set')
    from sklearn.datasets import make_moons
    Xmoon, Ymoon = make_moons(200, noise=0.05, random_state=0)
    plt.close()
    fig, ax = plt.subplots(1,1)
    
    x = Xmoon[:,0]
    y = Xmoon[:,1]
    ax.scatter(x,y)
    plt.show()

    max_components = 21
    n_components = np.arange(1,max_components)
    models = [GaussianMixture(n_components=n,covariance_type='full',random_state=0).fit(Xmoon) 
            for n in n_components]

    bic, idx = min((val, idx) for (idx, val) in enumerate([m.bic(Xmoon) for m in models]))
    bic_n_components = n_components[idx]

    aic, idx = min((val, idx) for (idx, val) in enumerate([m.aic(Xmoon) for m in models]))
    aic_n_components = n_components[idx]

    print("bic:{}@{} componenta".format(bic,bic_n_components))
    print("aic:{}@{} componenta".format(aic,aic_n_components))
    plt.close()
    plt.plot(n_components, [m.bic(Xmoon) for m in models], label='BIC')
    plt.plot(n_components, [m.aic(Xmoon) for m in models], label='AIC')
    plt.legend(loc='best')
    plt.xlabel('n_components')
    plt.show()


