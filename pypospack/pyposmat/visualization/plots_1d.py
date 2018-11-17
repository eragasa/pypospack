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
from scipy.stats import norm,gaussian_kde
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization

class PyposmatErrorNormalizationError(Exception): pass

class Pyposmat1DHistogramWithDensityPlots(PyposmatDataFileVisualization):

    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)

        self.plot_is_transparent = True
        self.plot_format = "eps"
   
    def read_configuration(self,filename):
        PyposmatDataFileVisualization.read_configuration(self,filename)

        self.qoi_targets = self.configuration.qoi_targets

    #def plot(self,x_name,results_filename_dict,filename=None):
    #    for k,v in results_filename_dict.items():
    #       print(k,v)

    def plot(self,
            x_name,
            x_min=None,
            x_max=None,
            x_step=None,  #DEPRECATED
            x_nsteps=1000,
            filename=None,
            include_histogram=True,
            include_normal=True,
            include_kde=True):

        self.histogram_nbins = 50
        self.histogram_color = "skyblue"
        self.histogram_label = "histogram"
        self.n_smallest_percentile = 0.95

        self.normal_dist_label = 'normal'
        self.normal_dist_color = 'black'
        self.normal_dist_linetype = '--'
        
        self.kde_label = 'KDE'
        self.kde_linetype = '-'
        self.kde_color = 'black'
        
        if x_name.endswith(".nerr"):
            self.calculate_normed_errors()
            
        if x_name in ["sum_all.nerr","sum_sq.nerr"]:
            _n_smallest_percentile = self.n_smallest_percentile
            _nrows,= self.df[x_name].shape
            _nsmallest = int(_nrows*_n_smallest_percentile)
            _x = self.df[x_name].nsmallest(_nsmallest)
            
        elif x_name == "sum_sq.nerr":
            _n_smallest_percentile = self.n_smallest_percentile
            _nrows,= self.df[x_name].shape
            _nsmallest = int(_nrows*_n_smallest_percentile)
            _x = self.df[x_name].nsmallest(_nsmallest)
        else:
            _x = self.df[x_name]

 
        if x_min is None:
            if x_name in ["sum_all.nerr","sum_sq.nerr"]: 
                _x_min = 0
            else: 
                _x_min = x.min()
        else: 
            _x_min = x_min
        
        if x_max is None: 
            _x_max = _x.max()
        else:
            _x_max = x_max
       

        self.fig, self.ax = plt.subplots()
        _fig = self.fig
        _ax = self.ax
        
        if include_histogram:
            # the histogram of the data
            _num_bins = self.histogram_nbins
            _histogram_color = self.histogram_color
            _histogram_label = self.histogram_label
            _n, _bins, _patches = _ax.hist(
                    _x, 
                    _num_bins, 
                    color = _histogram_color,
                    normed=1,
                    label= _histogram_label)
        
        if any([include_normal,include_kde]):
            _x_nsteps = x_nsteps
            _X = np.linspace(_x_min,_x_max,_x_nsteps)
        
        if include_normal: 

            _normal_dist_label = self.normal_dist_label
            _normal_dist_color =  self.normal_dist_color
            _normal_dist_linetype = self.normal_dist_linetype

            _mu = _x.mean(axis=0)
            _sigma = _x.std(axis=0)

            handle_norm = _ax.plot(
                    _X,
                    norm.pdf(X,_mu,_sigma),
                    _normal_dist_linetype,
                    label=_normal_dist_label,
                    color=_normal_dist_color)
        
        if include_kde:

            _kde_label = self.kde_label
            _kde_linetype = self.kde_linetype
            _kde_color = self.kde_color
            kde = gaussian_kde(_x)
            handle_kde = _ax.plot(
                    _X,
                    kde(_X),
                    _kde_linetype,
                    label=_kde_label,
                    color = _kde_color)
        
        # create_legend()
        if x_name == "sum_all.nerr":
            x_label = "Sum of Standardized Absolute Errors"
        elif x_name == "sum_sq.nerr":
            x_label = "Sum of Square Differences"
        else:
            x_label = self.configuration.latex_labels[x_name]
        y_label = "Probability Density"
       
        _ax.set_xlim(_x_min,_x_max)
        _ax.legend()
        _ax.set_xlabel(x_label)
        _ax.set_ylabel(y_label)

        if filename is not None:
            self.save_plot(fig=_fig,filename=filename)

    def save_plot(self,fig,filename):
        if fig is not None: 
            self.fig = fig
        _fig = fig
            
        # Tweak spacing to prevent clipping of ylabel
        _plot_is_transparent = self.plot_is_transparent
        _plot_format = self.plot_format
        _fig.tight_layout()
        _fig.show()
        _fig.savefig(
                filename,
                format = 'eps', 
                transparent=True)

        

    def calculate_normed_errors(self,df=None,qoi_names=None):
        """

        If a pandas.DataFrame is passed to df, then it will be set as the df attribute 
        for this class.  It willl
        Args:
            df (pandas.DataFrame)
            qoi_names (list) - a list of qoi names.  Default behavior will use all the 
                qois specified in the configuration object

        """
        if df is not None:
            self.df = copy.deepcopy(df)

        if qoi_names is not None:
            _qoi_names = list(qoi_names)
        else:
            _qoi_names = list(self.qoi_names)

        self.normed_error_names = []
        self.normed_error_validation_names = []
        for qn in _qoi_names:
            
            if qn in self.qoi_names:

                en = "{}.err".format(qn)
                nen = "{}.nerr".format(qn)
                self.normed_error_names.append(nen)
                
                
                q = self.qoi_targets[qn]

            elif qn in self.qoi_validation_names:
                en = "{}.err_v".format(qn)
                nen = "{}.nerr_v".format(qn)
                self.normed_error_validation_names.append(nen)

                q = self.qoi_validation_targets[qn]
            else:
                s = 80*"-"+"\n"
                s += "{:^80}\n".format('debugging information')
                s += 80*"-"+"\n"
                s += "qoi_name:{}\n".format(qn)
                s += "qoi_names\n"
                s += "\n".join(["  {}".format(v) for v in _qoi_names])+"\n"
                s += 80*"-"+"\n"
                s += "{:^80}\n".format('debugging information')
                s += "qoi_names\n"
                print(s)
                raise ValueError()

            self.df[nen] = self.df[qn]/q-1
        # sum normed errors
        self.df["sum_all.nerr"] = self.df[self.normed_error_names].abs().sum(axis=1)
        
        # TODO: do this in one-line
        _temp_df = self.df[self.normed_error_names]**2
        self.df["sum_sq.nerr"] = _temp_df.sum(axis=1)

        assert "sum_sq.nerr" in self.df
if __name__ == "__main__":
    from numpy.linalg.linalg import LinAlgError

    data_directory = "../../../data/MgO_pareto_data"
    fn_config = os.path.join(data_directory,'pyposmat.config.in')
    fn_results = os.path.join(data_directory,"results_000.out")

    myplot = Pyposmat1DHistogramWithDensityPlots()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    
    for pn in myplot.parameter_names:
        plot_fn = "{}.eps".format(pn)
        try:

            x_min = 1.2
            x_max = 2.7
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
