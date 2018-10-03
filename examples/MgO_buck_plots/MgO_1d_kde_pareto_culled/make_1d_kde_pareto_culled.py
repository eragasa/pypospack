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

import os,copy
from collections import OrderedDict
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats import norm,gaussian_kde
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization

class Pyposmat1DHistogramWithDensityPlots(PyposmatDataFileVisualization):

    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)

    def plot(self,x_name,results_filename_dict,filename=None):
        for k,v in results_filename_dict.items():
           print(k,v)

    def read_datafiles(self,results_fn,pareto_fn,culled_fn):

        self.results_datafile = PyposmatDataFile()
        self.results_datafile.read(results_fn)

        self.pareto_datafile = PyposmatDataFile()
        self.pareto_datafile.read(pareto_fn)

        self.culled_datafile = PyposmatDataFile()
        self.culled_datafile.read(culled_fn)

        self._parameter_names = self.results_datafile.parameter_names
        self._qoi_names = self.results_datafile.qoi_names
        self._error_names = self.results_datafile.error_names

        self.results_df = copy.deepcopy(self.results_datafile.df)
        self.pareto_df = copy.deepcopy(self.pareto_datafile.df)
        self.culled_df = copy.deepcopy(self.culled_datafile.df)

        for q in self.qoi_names:
            aen = '{}.abserr'.format(q)
            en = '{}.err'.format(q)
            self.results_df[aen] = self.results_df[en].abs()
            self.pareto_df[aen] = self.pareto_df[en].abs()
            self.culled_df[aen] = self.culled_df[en].abs()

    def plot(self,
            x_name,
            x_min=None,
            x_max=None,
            x_step=None,
            filename=None):

        fig, ax = plt.subplots()

        if x_min is None: 
            x_min = min(
                    self.results_df[x_name].min(),
                    self.pareto_df[x_name].min(),
                    self.culled_df[x_name].min()
                )
        if x_max is None: 
            x_max = max(
                    self.results_df[x_name].max(),
                    self.pareto_df[x_name].max(),
                    self.culled_df[x_name].max()
                )
        if x_step is None: 
            x_step = 0.01
        
        print("  x_range:{},{}".format(x_min,x_max))     
        X = np.arange(x_min,x_max,x_step)

        results_kde = gaussian_kde(self.results_df[x_name])
        results_hdl= ax.plot(
                X,
                results_kde(X),
                '--',
                label='results')

        pareto_kde = gaussian_kde(self.pareto_df[x_name])
        pareto_hdl = ax.plot(
                X,
                pareto_kde(X),
                '--',
                label='pareto')

        culled_kde = gaussian_kde(self.culled_df[x_name])
        pareto_hdl = ax.plot(
                X,
                culled_kde(X),
                '--',
                label='culled')

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
    from numpy.linalg.linalg import LinAlgError

    data_directory = "../../../data/MgO_pareto_data"
    fn_config = os.path.join(data_directory,'pyposmat.config.in')
    
    _results_fn= os.path.join(data_directory,"results_000.out")
    _pareto_fn = os.path.join(data_directory,"pareto_000.out")
    _culled_fn = os.path.join(data_directory,"culled_000.out")

    myplot = Pyposmat1DHistogramWithDensityPlots()
    myplot.read_datafiles(
            results_fn=_results_fn,
            pareto_fn=_pareto_fn,
            culled_fn=_culled_fn)
    myplot.read_configuration(filename=fn_config)
    
    for pn in myplot.parameter_names:
        plot_fn = "{}.eps".format(pn)
        try:

            x_min = 1.2
            x_max = 2.8
            myplot.plot(
                x_name=pn,
                x_min=x_min,
                x_max=x_max,
                filename=plot_fn
            )
            msg = " saving kde plot of, {}, to {}".format(pn,plot_fn)
            print(msg)
        except ZeroDivisionError as e:
            msg = "cannot plot the variable, {}, because the variable is deterministic".format(pn)
            print(msg)
        except LinAlgError as e:
            msg = "cannot plot the variable, {}, because the variable is deterministic".format(pn)
            print(msg)
