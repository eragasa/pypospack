import numpy as np
import scipy.stats

# References for KDE
# https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.gaussian_kde.html

_filename_png = 'kde_1d.png'
# caibrate the kde
mu = 0.0
sigma = 1.0
_loc = mu
_scale = sigma
_norm_rv = scipy.stats.norm(loc=_loc,scale=_scale)
data = _norm_rv.rvs(size=100)

print('simulation_data:{}'.format(str(type(data))))
print(data)

kde = scipy.stats.gaussian_kde(data)
print('kde_object:{}'.format(str(type(kde))))

# create the range of value for which we want to plot the kde
x_min = -5
x_max = 5
x_npoints = 1000
x = np.linspace(x_min,x_max,x_npoints)

# value of the KDE evaluate on the grid defined by x
kde_pdf_x = kde(x)
print('kde_pdf_x:{}'.format(str(type(kde_pdf_x))))

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,kde_pdf_x)
fig.savefig(_filename_png)

