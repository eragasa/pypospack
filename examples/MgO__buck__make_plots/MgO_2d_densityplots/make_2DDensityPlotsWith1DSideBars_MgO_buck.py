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
from scipy.stats import gaussian_kde

class Pyposmat2DDensityPlotsWith1DSidebars(PyposmatDataFileVisualization):
    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)
   
    def make_kde(self,x,y,h=None):

        values=np.vstack([x,y])
        if h is None:
            kde = gaussian_kde(values)
        elif h in ['silverman','silverman1986']:
            kde = gaussian_kde(values,_h)
        elif h is 'chiu1999':
            _h = Chiu1999_h(values)
            kde = gaussian_kde(values,_h)
        else:
            hde = gaussian_kde(values,h)
        return kde

    def set_axes_labels(self,ax,x_name,y_name,fontsize=25):
        ax.set_xlabel(x_name,fontsize=fontsize)
        ax.set_ylabel(y_name,fontsize=fontsize)

    def get_data_limits(self,name):
        xmin = self.df[name].min()
        xmax = self.df[name].max()
        return xmin,xmax
    
    def ax_scatter(self,
            ax,
            x : np.ndarray,
            y : np.ndarray,
            xy_marker_type='.',
            xy_marker_color='k',
            xy_marker_size=2,):

        ax.plot(
            x,
            y,
            '{}{}'.format(
                xy_marker_color,
                xy_marker_type),
            markersize=xy_marker_size)

    def ax_density_plot(self,
            ax,
            Z,
            XY_density_h,
            XY_cmap_name,
            xmin,
            xmax,
            ymin,
            ymax
        ):
        aspectratio = (xmax-xmin)/(ymax-ymin)
        ax.imshow(
            np.rot90(Z),
            cmap=plt.get_cmap(XY_cmap_name),
            extent=[xmin,xmax,ymin,ymax],
            aspect=aspectratio
            )
    
    def plot(self,
            x_name,
            y_name,
            x_lims = None,
            y_lims = None,
            show_XY_density=True,
            ax_name = None,
            XY_density_h = None,
            XY_cmap_name='Blues',
            fn_plot_out=None,
            show_xy_scatter=True,
            xy_marker_type='.',
            xy_marker_size=1,
            xy_marker_color='k',
            plot_fn = None):
        
        #plt.rc('text', usetex=True)

        # Define the x and y data
        x = self.df[x_name]
        y = self.df[y_name]

        if x_lims is None:
            xlims = [min(x),max(x)]
        else:
            xlims = xlims
        
        if y_lims is None:
            ylims = [min(y),max(y)]
        else:
            ylims = ylims
            

        # Define the locations for the axes
        left, width = 0.12, 0.55
        bottom, height = 0.12, 0.55
        bottom_h = bottom + height + 0.02
        left_h = left + width + 0.02
        # Set up the geometry of the three plots
        rect_2Dpane =   [left,   bottom,   width, height] # dimensions of temp plot
        rect_1Dpane_x = [left,   bottom_h, width, 0.25  ] # dimensions of x-histogram
        rect_1Dpane_y = [left_h, bottom,   0.25,  height] # dimensions of y-histogram

        # Set up the size of the figure
        fig = plt.figure(1, figsize=(9,9))

        # Make the three plots
        ax_2Dpane_XY = plt.axes(rect_2Dpane) # temperature plot
        ax_1Dpane_x = plt.axes(rect_1Dpane_x) # x histogram
        ax_1Dpane_y = plt.axes(rect_1Dpane_y) # y histogram

        # Remove the inner axes numbers of the histograms
        nullfmt = NullFormatter()
        ax_1Dpane_x.xaxis.set_major_formatter(nullfmt)
        ax_1Dpane_y.yaxis.set_major_formatter(nullfmt)

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
        kde = self.make_kde(x,y,h=XY_density_h)
        X,Y = np.mgrid[xmin:xmax:100j,ymin:ymax:100j]
        XY_grid = np.vstack([X.ravel(),Y.ravel()])
        Z=np.reshape(kde(XY_grid),X.shape)

        if show_XY_density is True:
            self.ax_density_plot(
                ax_2Dpane_XY,
                Z,
                XY_density_h,
                XY_cmap_name,
                xmin,
                xmax,
                ymin,
                ymax)

        if show_xy_scatter is True:
            self.ax_scatter(
                ax_2Dpane_XY,
                x,
                y,
                xy_marker_type=xy_marker_type,
                xy_marker_color=xy_marker_color,
                xy_marker_size=xy_marker_size)

        #Plot the axes labels
        xlabel = self.configuration.latex_labels[x_name]
        ylabel = self.configuration.latex_labels[y_name]
        
        print('xlabel:{}'.format(xlabel))
        print('ylabel:{}'.format(ylabel))
        
        ax_2Dpane_XY.set_xlabel(xlabel,fontsize=12)
        ax_2Dpane_XY.set_ylabel(ylabel,fontsize=12)

        #Make the tickmarks pretty
        ticklabels = ax_2Dpane_XY.get_xticklabels()
        for label in ticklabels:
            label.set_fontsize(8)
            label.set_family('serif')

        ticklabels = ax_2Dpane_XY.get_yticklabels()
        for label in ticklabels:
            label.set_fontsize(8)
            label.set_family('serif')

        #Set up the plot limits
        ax_2Dpane_XY.set_xlim(xlims)
        ax_2Dpane_XY.set_ylim(ylims)

        #Set up the histogram bins
        xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
        ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)

        #Plot the histograms
        is_plot_1D_histograms = False
        if is_plot_1D_histograms:
            ax_1Dpane_x.hist(x, bins=xbins, color = 'blue')
            ax_1Dpane_y.hist(y, bins=ybins, orientation='horizontal', color = 'red')

        #Plot the KDEs
        kde_x = gaussian_kde(x)
        kde_y = gaussian_kde(y)
        _x = np.linspace(min(x),max(x),1000)
        _y = np.linspace(min(y),max(y),1000)
        ax_1Dpane_x.plot(_x,kde_x(_x))
        #ax_1Dpane_x.plot(_x,kde_x(_x))
        ax_1Dpane_y.plot(kde_y(_y),_y)
 
        #Set up the panel limits limits
        ax_1Dpane_x.set_xlim( min(x), max(x) )
        ax_1Dpane_y.set_ylim( min(y), max(y) )

        #Make the tickmarks pretty
        ticklabels = ax_1Dpane_x.get_yticklabels()
        for label in ticklabels:
            label.set_fontsize(8)
            label.set_family('serif')

        #Make the tickmarks pretty
        ticklabels = ax_1Dpane_y.get_xticklabels()
        for label in ticklabels:
            label.set_fontsize(8)
            label.set_family('serif')

        #Cool trick that changes the number of tickmarks for the histogram axes
        ax_1Dpane_y.xaxis.set_major_locator(MaxNLocator(4))
        ax_1Dpane_x.yaxis.set_major_locator(MaxNLocator(4))

        # Save to a File
        if plot_fn is None:
            plot_fn = '{}__{}.eps'.format(x_name,y_name)
        msg = "saving plot as {}".format(plot_fn)
        plt.savefig(plot_fn, transparent=True)

if __name__ == "__main__":
    data_directory = "../../../data/MgO_pareto_data"
    fn_config = os.path.join(data_directory,'pyposmat.config.in')
    fn_results = os.path.join(data_directory,"culled_009.out")

    
    x_name = 'OO_A'
    y_name = 'OO_rho'
    myplot = Pyposmat2DDensityPlotsWith1DSidebars()
    myplot.read_datafile(filename=fn_results)
    myplot.read_configuration(filename=fn_config)
    myplot.plot(x_name=x_name,y_name=y_name)
