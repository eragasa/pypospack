import copy,os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class PyposmatDatafileVisualization(object):
    def __init__(self):
        self.datafile = None
        self.configuation = None
    @property
    def parameter_names(self):
        return self._parameter_names

    @property
    def qoi_names(self):
        return self._qoi_names

    @property
    def error_names(self):
        return self._error_names

    @property
    def df(self):
        return self._df

    @property
    def parameter_df(self):
        return self._parameter_df

    @property
    def qoi_df(self):
        return self._qoi_df

    @property
    def error_df(self):
        error_df(self)
        return self._error_df

    def read_datafile(self,filename):
        self.datafile = PyposmatDataFile()
        self.datafile.read(filename=filename)

        self._parameter_names = self.datafile.parameter_names
        self._qoi_names = self.datafile.qoi_names
        self._error_names = self.datafile.error_names

        self._df = copy.deepcopy(self.datafile.df)
        self.create_absolute_errors()

    def create_absolute_errors(self):
        for q in self.qoi_names:
            aen = "{}.abserr".format(q)
            en = "{}.err".format(q)
            self._df[aen] = self._df[en].abs()

    def read_configuration(self,filename):
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)

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

    def set_labels(self,ax,x_name,y_name):
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
class Pyposmat3DScatterWithProjections(PyposmatDatafileVisualization):

    def __init__(self):
        PyposmatDatafileVisualization.__init__(self)

    def plot(self,x_name,y_name,z_name,filename=None):
        x= self.df[x_name].abs()
        y= self.df[y_name].abs()
        z= self.df[z_name].abs()

        x_min = x.min()
        x_max = x.max()
        y_min = y.min()
        y_max = y.max()
        z_min = z.min()
        z_max = z.max()

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        ax.plot(x, z, 'r+', zdir='y', marker='.',markersize=1)
        ax.plot(y, z, 'g+', zdir='x', marker='.',markersize=1)
        ax.plot(x, y, 'k+', zdir='z', marker='.',markersize=1)
        ax.scatter(x,y,z,marker='.',s=1,c='k')

        self.set_labels(ax,x_name,y_name,z_name)
        #ax.set_xlim([x_min,x_max])#ax.set_xlim([-0.5, 1.5])
        #ax.set_ylim([y_min,y_max])#ax.set_ylim([-0.5, 1.5])
        #ax.set_zlim([z_min,z_max])#ax.set_zlim([-1.5, 1.5])

        if filename is None:
            plt.show()
        else:
            fig.savefig(filename)

    def set_labels(self,ax,x_name,y_name,z_name):
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_zlabel(z_name)
if __name__ == "__main__":
    fn_config = "pyposmat.config.in"
    fn_results = "pyposmat.kde.out"
    myplot = Pyposmat3DScatterWithProjections()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    myplot.plot(
        x_name='Ni_fcc.c11.err',
        y_name='Ni_fcc.c12.err',
        z_name='Ni_fcc.c44.err'
    )

    myplot = Pyposmat2DScatterPlot()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    myplot.is_pareto = True
    myplot.plot(
        x_name='Ni_fcc.c11.abserr',
        y_name='Ni_fcc.c12.abserr'
    )
