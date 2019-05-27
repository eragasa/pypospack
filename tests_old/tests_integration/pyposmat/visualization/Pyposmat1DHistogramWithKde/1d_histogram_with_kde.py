"""
=========================================================
Demo of the histogram (hist) function with a few features
=========================================================

In addition to the basic histogram, this demo shows a few optional
features:

    * Setting the number of data bins
    * The ``normed`` flag, which normalizes bin heights so that the
      integral of the histogram is 1. The resulting histogram is an
      approximation of the probability density function.
    * Setting the face color of the bars
    * Setting the opacity (alpha value).

Selecting different bin counts and sizes can significantly affect the
shape of a histogram. The Astropy docs have a great section on how to
select these parameters:
http://docs.astropy.org/en/stable/visualization/histogram.html
"""

import os
from collections import OrderedDict
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization

class Pyposmat1DHistogramWithDensityPlots(PyposmatDataFileVisualization):

    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)

    def plot(self,x_name,filename=None):
        x = self.df[x_name]
        mu = x.mean(axis=0)
        sigma = x.std(axis=0)

        num_bins = 50

        fig, ax = plt.subplots()

        # the histogram of the data
        n, bins, patches = ax.hist(
                x, 
                num_bins, 
                color = 'skyblue',
                normed=1,
                label='histogram')

        # add a 'best fit' line
        y = mlab.normpdf(bins, mu, sigma)
        
        kde = gaussian_kde(x)
        X = np.arange(
                self.df[x_name].min(),
                self.df[x_name].max(),0.01)

        handle_norm = ax.plot(bins, y, '--',label='normal')
        handle_kde = ax.plot(X,kde(X),'-',label='KDE',color = 'black')
        ax.legend()
        ax.set_xlabel(r'$\beta$')
        ax.set_ylabel('Probability density')

        # Tweak spacing to prevent clipping of ylabel
        fig.tight_layout()
        fig.show()
        fig.savefig(filename,format = 'eps', transparent=True)


if __name__ == "__main__":
    data_dir = "../../../../data_test/Ni__eam__born_exp_fs_00/data__Ni__eam__born_exp_fs_02" 
    fn_config = os.path.join(data_dir,"pyposmat.config.in")
    fn_results = os.path.join(data_dir,"pyposmat.kde.9.out")
    myplot = Pyposmat1DHistogramWithDensityPlots()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    for pn in myplot.parameter_names:
        myplot.plot(
            x_name=pn,
            filename=pn+".eps"
        )
