import numpy as np
import scipy.stats

# testing for normal distribution

from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

from pypospack.statistics import kullbach_lieber_divergence
def generate_spd_matrix(n):
    """
    Generate an nxn symmetric, positive definite matrix

    Args:

    n(int): 
    """
    A = np.random.rand(n,n)

    A = 0.5*(A+A.T)

    A = A * n*np.eye(n)

    return A

n_samples_normal = 1000
n_samples_kde = 1000

rv_norm = norm(0,1)
X_norm = rv_norm.rvs(size=1000)
print('X_norm',X_norm.shape)

rv_kde_1 = gaussian_kde(X_norm)
X_kde = rv_kde_1.resample(size=1000)
print('X_kde',X_kde.shape)

rv_kde_2 = gaussian_kde(X_kde)
kld = kullbach_lieber_divergence(rv_kde_1,rv_kde_2,1000)

print(kld)
exit()
import matplotlib.pyplot as plt
xmin = min(X_norm.min(),X_kde.min())
xmax = max(X_norm.max(),X_kde.max())
x = np.linspace(xmin,xmax,1000)
fig, ax = plt.subplots()
ax.plot(x,rv_norm.pdf(x))
ax.plot(x,rv_kde.pdf(x))
plt.show()


