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
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import gaussian_kde
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization

class Pyposmat1DKdeDensityPlots(PyposmatDataFileVisualization):

    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)
        self.n_iterations = None

    def get_linestyle(self,i=None,name=None):

        assert i is None or isinstance(i,int)
        assert name is None or isinstance(name,str)


        linestyles = OrderedDict([
            ('solid',               (0, ())),
            #('loosely dotted',      (0, (1, 10))),
            #('dotted',              (0, (1, 5))),
            ('densely dotted',      (0, (1, 1))),
            
            #('loosely dashed',      (0, (5, 10))),
            #('dashed',              (0, (5, 5))),
            ('densely dashed',      (0, (5, 1))),
            
            #('loosely dashdotted',  (0, (3, 10, 1, 10))),
            #('dashdotted',          (0, (3, 5, 1, 5))),
            ('densely dashdotted',  (0, (3, 1, 1, 1))),
            
            #('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
            #('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
            ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]
        )

        if isinstance(i,int):
            for i_linestyle, (k,v) in enumerate(linestyles.items()):
                if i_linestyle == i:
                    return v
        elif isinstance(name,str):
            for k,v in linestyles.items():
                if k == name:
                    return v
           
    def plot(self,
            x_name, 
            x_min=None,
            x_max=None,
            x_step=None,
            filename=None,
            plt_labels=None,
            colormap='cool',
            show_iterations=None):

        """

            show_iterations (list of int)
        """
        fig, ax = plt.subplots()

        plot_handles = OrderedDict()

        if type(self.df) is OrderedDict:

            if x_min is None:
                x_min = min([df[x_name].min() for k,df in self.df.items()])
            if x_max is None:
                x_max = max([df[x_name].max() for k,df in self.df.items()])
            if x_step is None:
                x_step = 0.01
      
            if plt_labels is None:
                if show_iterations is not None:
                    plt_labels = ['N={}'.format(i) for i in range(self.n_iterations)]
                else:
                    plt_labels = ['N={}'.format(i) for i in show_iterations]
            
            if type(colormap) is str:
                _cm_cmap = plt.get_cmap(colormap)
            elif type(colormap) is list:

                # make custom gradient colormap
                import matplotlib.colors as mcol
                _cm_cmap_name = "{}To{}".format(colormap[0],colormap[1])
                _cm_cmap = mcol.LinearSegmentedColormap.from_list(_cm_cmap_name,colormap)

            else:
                raise ValueError()
            _cm_norm = colors.Normalize(vmin=0,vmax=len(self.df))
            _cm_map = cm.ScalarMappable(
                    norm = _cm_norm,
                    cmap = _cm_cmap)


            if show_iterations is None:
                _show_iterations = [k for k in self.df]
            else:
                _show_iterations = show_iterations

            print(_show_iterations)
            i_linestyle = 0
            for i, (k,df) in enumerate(self.df.items()):
                if i in _show_iterations:
                    print('iteration {}: working from file {}, n_rows:{}'.format(i,k,df[x_name].shape[0]))

                    X = np.arange(x_min,x_max,x_step)
                    kde = gaussian_kde(df[x_name])
                    linestyle = self.get_linestyle(i_linestyle)
                    color = _cm_map.to_rgba(i)

                    plot_handles[k] = ax.plot(
                            X,
                            kde(X),
                            linestyle=self.get_linestyle(i_linestyle),
                            label=plt_labels[i],
                            color = _cm_map.to_rgba(i)
                    )
                    i_linestyle += 1
                else:
                    print('iteration {}: skipping'.format(i))

        elif type(self.df) is pd.DataFrame:
            print(self.filename)
        else:
            errmsg("neither the filename nor the filenames attribute has been set")
            raise ValueError(errmsg)

        # create_legend()
        x_label = self.configuration.latex_labels[x_name]
        y_label = "Probability Density"
        
        ax.legend()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        ax.set_xlim(x_min,x_max)  
        
        # Tweak spacing to prevent clipping of ylabel
        fig.tight_layout()
        fig.show()
        fig.savefig(filename,format = 'eps', transparent=True)

def get_result_filenames(data_dir,fn_format,n_iterations):
    fn_results = []
    for i in range(n_iterations):
        fn_results.append(
                os.path.join(data_directory,fn_format.format(i))
            )

    return fn_results

def validate_result_filenames(fn_results):
    if not isinstance(fn_results,list):
        err_msg = "type(fn_results):{}\n".format(type(fn_results))
        err_msg += "\tshould be {}".format(type(list))
        raise ValueError(err_msg)
    for i,fn in enumerate(fn_results):
        print("{}:{}".format(i,fn))

if __name__ == "__main__":
    data_directory = "../../../data/MgO_pareto_data"
    fn_config = os.path.join(data_directory,'pyposmat.config.in')
    result_fn_format = "results_{:03d}.out"
    n_iterations = 10
    colormap=['red','blue']


    show_iterations=[0,1,2,4,9]
    fn_results = get_result_filenames(
            data_dir=data_directory,
            fn_format=result_fn_format,
            n_iterations=n_iterations
        )

    validate_result_filenames(fn_results)

    myplot = Pyposmat1DKdeDensityPlots()
    myplot.read_datafiles(filenames=fn_results)
    myplot.read_configuration(filename=fn_config)
   
    from numpy.linalg.linalg import LinAlgError
    for pn in myplot.parameter_names:
        print("trying to create the plot for parameter {}".format(pn))
        plot_fn = "{}.eps".format(pn)
        try:
            myplot.plot(x_name=pn,filename=plot_fn,colormap = colormap, show_iterations=show_iterations)
            msg = "saving kde plot of, {}, to {}".format(pn,plot_fn)
            print(msg)
        except ZeroDivisionError as e:
            msg = "cannot plot the variable, {}, because the variable is deterministic".format(pn)
            print(msg)
        except LinAlgError as e:
            msg = "cannot plot the variable, {}, because the variable is deterministic".format(pn)
            print(msg)
