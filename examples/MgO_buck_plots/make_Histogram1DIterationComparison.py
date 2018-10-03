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
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization

class Pyposmat1DHistogramWithDensityPlots(PyposmatDataFileVisualization):

    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)

    def plot_iterative_comparison(self,x_name,results_filename_dict,filename=None):
       
        msg = "reading in the datafiles..."
        print(msg)
        self.datafiles = OrderedDict()
        for k,v in results_filename_dict.items():
           print(k,v)
           self.datafiles[k] = PyposmatDataFile()
           self.datafiles[k].read(v)

        for i in results_filename_dict:
            x = self.datafiles[i].df[x_name]
            kde = gaussian_kde(x)

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
        
        # create_legend()
        x_label = self.configuration.latex_labels[x_name]
        y_label = "Probability Density"
        
        ax.legend()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        # Tweak spacing to prevent clipping of ylabel
        fig.tight_layout()
        fig.show()
        fig.savefig(filename,format = 'eps', transparent=True)


if __name__ == "__main__":
    data_directory = "../../data/MgO_pareto_data"
    fn_config = os.path.join(data_directory,'pyposmat.config.in')
    fn_results = os.path.join(data_directory,"culled_005.out")

    n_iterations = 5
    fn_results_dict = OrderedDict()
    for i in range(n_iterations):
        fn_results_dict[i] = os.path.join(
                data_directory,
                "culled_00{}.out".format(i))

    myplot = Pyposmat1DHistogramWithDensityPlots()
    myplot.read_configuration(filename=fn_config)
    myplot.plot_iterative_comparison(
            x_name='chrg_Mg',
            results_filename_dict=fn_results_dict,
            filename='chrg_Mg.eps')
    #for pn in myplot.parameter_names:
    #    plot_fn = "{}.eps".format(pn)
    #    try:
    #        myplot.plot(
    #            x_name=pn,
    #            filename=plot_fn
    #        )
    #        msg = " saving kde plot of, {}, to {}".format(pn,plot_fn)
    #        print(msg)
    #    except ZeroDivisionError as e:
    #        msg = "cannot plot the variable, {}, because the variable is deterministic".format(pn)
    #        print(msg)
