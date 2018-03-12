import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace

import copy,os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
from matplotlib import rc

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization

class Pyposmat2DDensityPlots(PyposmatDataFileVisualization):
    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)
    
    def plot(self,x_name,y_name,filename=None,
            xy_cmap_name='Blues'):
        
        plt.rc('text', usetex=True)
        #plt.ion()


        # Define a function to make the ellipses

        # Define the x and y data
        # For example just using random numbers
        x = self.df[x_name]
        y = self.df[y_name]

        # Set up default x and y limits
        xlims = [min(x),max(x)]
        ylims = [min(y),max(y)]

        # Set up your x and y labels
        xlabel = r'$a_{0,\phi} \mathrm{[\AA]}$'
        ylabel = r'$a_{a,\rho} \mathrm{[\AA]}$'

        # Define the locations for the axes
        left, width = 0.12, 0.55
        bottom, height = 0.12, 0.55
        bottom_h = left_h = left+width+0.02

        # Set up the geometry of the three plots
        rect_temperature = [left,   bottom,   width, height] # dimensions of temp plot
        rect_histx =       [left,   bottom_h, width, 0.25  ] # dimensions of x-histogram
        rect_histy =       [left_h, bottom,   0.25,  height] # dimensions of y-histogram

        # Set up the size of the figure
        fig = plt.figure(1, figsize=(9.5,9))

        # Make the three plots
        axHeatplot_xy = plt.axes(rect_temperature) # temperature plot
        axHistx = plt.axes(rect_histx) # x histogram
        axHisty = plt.axes(rect_histy) # y histogram

        # Remove the inner axes numbers of the histograms
        nullfmt = NullFormatter()
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        # Find the min/max of the data
        xmin = min(xlims)
        xmax = max(xlims)
        ymin = min(ylims)
        ymax = max(ylims)

        # Make the 'main' temperature plot
        # Define the number of bins
        nxbins = 50
        nybins = 50
        nbins = 100

        xbins = linspace(start = xmin, stop = xmax, num = nxbins)
        ybins = linspace(start = ymin, stop = ymax, num = nybins)
        xcenter = (xbins[0:-1]+xbins[1:])/2.0
        ycenter = (ybins[0:-1]+ybins[1:])/2.0
        aspectratio = 1.0*(xmax - xmin)/(1.0*ymax - ymin)

        Z, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins))
        X = xcenter
        Y = ycenter

        # Plot the temperature data
        axHeatplot_xy.imshow(
                Z,
                extent=[xmin,xmax,ymin,ymax],
                interpolation='nearest', 
                cmap=plt.get_cmap(xy_cmap_name),
                origin='lower',
                aspect=aspectratio)
         
        #Plot the axes labels
        axHeatplot_xy.set_xlabel(xlabel,fontsize=25)
        axHeatplot_xy.set_ylabel(ylabel,fontsize=25)

        #Make the tickmarks pretty
        ticklabels = axHeatplot_xy.get_xticklabels()
        for label in ticklabels:
            label.set_fontsize(18)
            label.set_family('serif')

        ticklabels = axHeatplot_xy.get_yticklabels()
        for label in ticklabels:
            label.set_fontsize(18)
            label.set_family('serif')

        #Set up the plot limits
        axHeatplot_xy.set_xlim(xlims)
        axHeatplot_xy.set_ylim(ylims)

        #Set up the histogram bins
        xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
        ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)

        #Plot the histograms
        axHistx.hist(x, bins=xbins, color = 'blue')
        axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'red')

        #Set up the histogram limits
        axHistx.set_xlim( min(x), max(x) )
        axHisty.set_ylim( min(y), max(y) )

        #Make the tickmarks pretty
        ticklabels = axHistx.get_yticklabels()
        for label in ticklabels:
            label.set_fontsize(12)
            label.set_family('serif')

        #Make the tickmarks pretty
        ticklabels = axHisty.get_xticklabels()
        for label in ticklabels:
            label.set_fontsize(12)
            label.set_family('serif')

        #Cool trick that changes the number of tickmarks for the histogram axes
        axHisty.xaxis.set_major_locator(MaxNLocator(4))
        axHistx.yaxis.set_major_locator(MaxNLocator(4))

        #Show the plot
        plt.draw()

        # Save to a File
        filename = 'myplot'
        plt.savefig(filename + '.eps',format = 'eps', transparent=True)

if __name__ == "__main__":
    data_dir = "../../../../data_test/Ni__eam__born_exp_fs_00/data__Ni__eam__born_exp_fs_02" 
    fn_config = os.path.join(data_dir,"pyposmat.config.in")
    fn_results = os.path.join(data_dir,"pyposmat.kde.9.out")
    myplot = Pyposmat2DDensityPlots()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    myplot.plot(
        x_name='p_NiNi_r0',
        y_name='d_Ni_r0'
    )
