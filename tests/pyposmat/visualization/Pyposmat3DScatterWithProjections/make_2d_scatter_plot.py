import copy,os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatDatafileVisualization

class Pyposmat2DScatterPlot(PyposmatDatafileVisualization):

    def __init__(self):
        PyposmatDatafileVisualization.__init__(self)

        self.is_pareto = False

    def plot(self,x_name,y_name,filename=None):
        fig = plt.figure(figsize=plt.figaspect(1))
        ax= fig.add_subplot(1,1,1)

        x = self.df[x_name]
        y = self.df[y_name]

        ax.scatter(x,y)

        self.set_labels(ax,x_name,y_name)

        self.set_axes_limits(ax)

        if filename is None:
            plt.show()
        else:
            fig.savefig(filename)

    def set_axes_limits(self,ax):

        if self.is_pareto:
            ax.set_xlim(left=0)
            ax.set_ylim(left=0)

    def get_latex_label(self,name):
        return name

    def set_labels(self,ax,x_name,y_name):
        x_latex_label = self.get_latex_label(x_name)
        y_latex_label = self.get_latex_label(y_name)

        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)


if __name__ == "__main__":
    fn_config = "pyposmat.config.in"
    fn_results = "pyposmat.kde.out"
    myplot = Pyposmat2DScatterPlot()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    myplot.is_pareto = True
    myplot.plot(
        x_name='Ni_fcc.c11.abserr',
        y_name='Ni_fcc.c12.abserr'
    )
