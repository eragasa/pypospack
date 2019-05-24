import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.stats import norm
#from pypospack.pyposmat.data import PyposmatConfigurationFile
#from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.visualization import PyposmatAbstractPlot

class PyposmatQoiPlot(PyposmatAbstractPlot):

    def __init__(self,
            config=None,
            data=None,
            is_plot_histogram = False,
            is_plot_kde = True):

        PyposmatAbstractPlot.__init__(self,config=config,data=data)

        self.qoi_name = None
        self.create_subplots()

        self.is_plot_histogram = False
        self.is_plot_kde = True
        self.is_plot_qoitarget = False

    def create_subplots(self):
        PyposmatAbstractPlot.create_subplots(self,nrows=1,ncols=1)

    def add_historgram(self,qoi_name):
        self.ax.hist(self.data.df[qoi_name])

    def add_kde_plot(self,qoi_name,label=None,x_limits=None,color=None):

        assert isinstance(qoi_name,str)

        if x_limits is None:
            x_lim_min = self.get_x_lim_min(x=self.data.df[qoi_name])
            x_lim_max = self.get_x_lim_max(x=self.data.df[qoi_name])
        else:
            x_lim_min = x_limits[0]
            x_lim_max = x_limits[1]

        x = np.linspace(x_lim_min,x_lim_max,1000)
        kde = gaussian_kde(self.data.df[qoi_name])

        if color is None:
            self.ax.plot(x,kde(x))
        else:
            self.ax.plot(x,kde(x),label=label,color=color)

    def add_qoitarget(self,qoi_name):
        self.ax.axvline(self.configuration.qoi_targets[qoi_name])

    def get_x_lim_min(self,x):

        if isinstance(x,pd.Series) or isinstance(x,np.ndarray):
            norm 
            return x.min()
        elif isinstance(x,list):
            return min(x)
        else:
            m = "type(x)={}".format(str(type(x)))
            raise TypeError(m)

    def get_x_lim_max(self,x):

        if isinstance(x,pd.Series) or isinstance(x,np.ndarray):
            return x.max()
        elif isinstance(x,list):
            return max(x)
        else:
            m = "type(x)={}".format(str(type(x)))
            raise TypeError(m)

    def add_qoi_plot(self,qoi_name,x_limits=None,label=None,color=None):

        if self.is_plot_histogram:
            self.add_histogram_plot(df)

        if self.is_plot_kde:
            self.add_kde_plot(qoi_name=qoi_name,x_limits=x_limits,label=label,color=color)

        if self.is_plot_qoitarget:
            self.plot_qoitarget(qoi_name_qoi_name)

    def legend(self,loc='upper right'):
        self.ax.legend(loc=loc)
