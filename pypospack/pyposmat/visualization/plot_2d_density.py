import numpy as np
import pandas as pd
from scipy import stats

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatDataFileVisualization  

class Pyposmat2DDensityPlots(PyposmatDataFileVisualization):
    def __init__(self):
        PyposmatDataFileVisualization.__init__(self)

        self.figure = None
        self.axes = None

    def make_kde(self,x,y,h=None):
        values=np.vstack([x,y])
        if h is None:
            kde = stats.gaussian_kde(values)
        elif h in ['silverman','silverman1986']:
            kde = stats.gaussian_kde(values,'silverman')
        elif h is 'chiu1999':
            _h = Chiu1999_h(values)
            kde = stats.gaussian_kde(values,_h)
        else:
            hde = stats.gaussian_kde(values,h)
        return kde

    def set_axes_labels(self,ax,x_name,y_name):
        
        try:
            x_label = self.configuration.latex_labels[x_name]['name']
        except KeyError as e:
            if e.args[0] == 'latex_labels':
                msg = 'the x_label is using {} because latex_label is not configured in the configuration file'.format(x_name)
                print(msg)
                x_label = x_name
            else:
                raise
        try:
            y_label = self.configuration.latex_labels[y_name]['name']
        except KeyError as e:
            if e.args[0] == 'latex_labels':
                msg = 'the y_label is using {} because latex_label is not configured in the configuration file'.format(y_name)
                print(msg)
                y_label = y_name
            else:
                raise

        print('x_label:{}'.format(x_label))
        print('y_label:{}'.format(y_label))
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

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
            x_name : str,
            y_name : str,
            x_lims = None, # list of length 2, x_lims[0] is the xmin, x_lims[1] is the xmax
            y_lims = None, # list of length 2, y_lims[0] is the ymin, y_lims[1] is the ymax
            fn_plot_out=True, # accepts either True or str
            ax_name = None,
            show_XY_density = True,
            XY_density_h = None,
            XY_cmap_name='Blues',
            show_xy_scatter=True,
            xy_marker_type='.',
            xy_marker_size=2,
            xy_marker_color='k',
            show_plot=False):
        
        x = self.df[x_name].values  # type: np.ndarray
        y = self.df[y_name].values  # type: np.ndarray
        xy_grid = np.vstack([x,y])
        kde = self.make_kde(x,y,h=XY_density_h)

        if x_lims is None:
            xmin,xmax = self.get_data_limits(x_name)
        else:
            xmin,xmax = x_lims
        
        if y_lims is None:
            ymin,ymax = self.get_data_limits(y_name)
        else:
            ymin,ymax = y_lims

        # create figure object if it doesn't exist
        if self.figure is None:
            self.figure = plt.figure()
        _fig = self.figure

        # creates axes dictionary if it doesn't exist
        if self.axes is None:
            self.axes = OrderedDict()

        # create new axes
        if ax_name is None:
            an = '{}__{}'.format(x_name,y_name)
        else:
            an = ax_name
        self.axes[an] = _fig.add_subplot(111)

        # pointer to actual object
        _ax = self.axes[an]

        # make surface
        X_density=100j
        Y_density=100j

        # build the grid
        X,Y = np.mgrid[xmin:xmax:100j,ymin:ymax:100j]
        XY_grid = np.vstack([X.ravel(),Y.ravel()])
        # evaluate density on the grid
        Z = np.reshape(kde(XY_grid),X.shape)
        
        if show_XY_density is True:
            #_ax.imshow(
            #    np.rot90(Z),
            #    cmap=plt.get_cmap(XY_cmap_name),
            #    extent=[xmin,xmax,ymin,ymax])
            self.ax_density_plot(
                _ax,
                Z,
                XY_density_h,
                XY_cmap_name,
                xmin,
                xmax,
                ymin,
                ymax)

        if show_xy_scatter is True:
            self.ax_scatter(
                    _ax,
                    x,
                    y,
                    xy_marker_type=xy_marker_type,
                    xy_marker_color=xy_marker_color,
                    xy_marker_size=xy_marker_size)

        self.set_axes_labels(_ax,x_name,y_name)
          
        if show_plot is True:
            plt.show()

        if fn_plot_out is not None:
            print("...writing out figure")
            if type(fn_plot_out) is str:
                fn = fn_plot_out
            elif fn_plot_out is True:
                fn = "{}__{}.png".format(x_name,y_name)
            elif fn_plot_out is False:
                pass
            print("...writing out figure_name:{}".format(fn))
            _fig.savefig(fn)
