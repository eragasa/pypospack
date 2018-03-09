import copy,os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatDatafileVisualization

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
