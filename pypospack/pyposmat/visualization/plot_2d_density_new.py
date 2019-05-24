import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pypospack.pyposmat.visualization import PyposmatAbstractPlot

class Pyposmat2DDensityPlot(PyposmatAbstractPlot):

    kde_bandwidth_types = ['silverman','silverman1986','chiu1999']

    def __init__(self,config=None,data=None):
        PyposmatAbstractPlot.__init__(self,config=config,data=data)

        self.x_limits = None
        self.y_limits = None
    def determine_limits(self,name,ppf_min=0.1,ppf_max=0.9):

        assert name in self.configuration.qoi_names \
               or name in self.configuration.parameter_names
        assert isinstance(ppf_min,float)
        assert isinstance(ppf_max,float)

        norm_rv =  stats.norm(
                loc = self.data.df[name].mean(),
                scale = self.data.df[name].std()
                )
        lim_min = norm_rv.ppf(ppf_min)
        lim_max = norm_rv.ppf(ppf_max)

        return lim_min,lim_max
    
    def determine_x_limits(self,x_name=None,x_limits=None,ppf_min=0.1,ppf_max=0.9):

        assert x_name is None or isinstance(x_name,str)
        assert x_limits is None or isinstance(x_name,list)
        assert isinstance(ppf_min,float)
        assert isinstance(ppf_max,float)

        if x_name is None:
            x_name = self.x_name

        if x_limits is None:
            x_lim_min,x_lim_max = self.determine_limits(x_name)
            self.x_limits = (x_lim_min,x_lim_max)
        else:
            self.x_limits = x_limits

        return self.x_limits

    def determine_y_limits(self,y_name=None,y_limits=None,ppf_min=0.1,ppf_max=0.9):

        assert y_name is None or isinstance(y_name,str)
        assert y_limits is None or isinstance(y_name,list)
        assert isinstance(ppf_min,float)
        assert isinstance(ppf_max,float)

        if y_name is None:
            y_name = self.y_name

        if y_limits is None:
            y_lim_min,y_lim_max = self.determine_limits(y_name)
            self.y_limits = (y_lim_min,y_lim_max)
        else:
            self.y_limits = y_limits

        return self.y_limits

    def plot(self,
             x_name,y_name,
             with_kde_plot=True,
             with_data_plot=True,
             x_limits=None,y_limits=None,h=None):

        assert x_name in self.configuration.qoi_names \
                or x_name in self.configuration.parameter_names
        assert y_name in self.configuration.qoi_names \
                or y_name in self.configuration.parameter_names
        assert x_limits is None \
                or isinstance(x_limits,list)
        assert y_limits is None \
                or isinstance(y_limits,list)
        assert h is None \
                or h in kde_bandwidth_types

        self.x_name = x_name
        self.y_name = y_name

        self.determine_x_limits()
        self.determine_y_limits()

        x = self.data.df[x_name].values
        y = self.data.df[y_name].values
      
        if self.fig is None or self.ax is None:
            self.create_subplots()

        if with_kde_plot:
            self.plot_kde(x,y,h)

        if with_data_plot:
            self.plot_data_points(x,y)

        self.ax.set_xlim(self.x_limits[0],self.x_limits[1])
        self.ax.set_ylim(self.y_limits[0],self.y_limits[1])

        xy_grid = np.vstack([x,y])
        kde = self.make_kde(x,y,h=h)

    def plot_kde(self,x,y,h=None,XY_cmap_name='Blues'):

        # build the grid
        xmin = self.x_limits[0]
        xmax = self.x_limits[1]
        ymin = self.y_limits[0]
        ymax = self.y_limits[1]
        X_density=200j
        Y_density=200j
        X,Y = np.mgrid[xmin:xmax:X_density,ymin:ymax:Y_density]
        XY_grid = np.vstack([X.ravel(),Y.ravel()])

        # evaluate density on the grid
        kde = self.make_kde(x,y,h)
        Z = np.reshape(kde(XY_grid),X.shape)

        aspectratio=(xmax-xmin)/(ymax-ymin)
        self.ax.imshow(
                np.rot90(Z),
                cmap=plt.get_cmap(XY_cmap_name),
                extent=[xmin,xmax,ymin,ymax],
                aspect=aspectratio)
    
    def plot_data_points(self,x,y,size=1):

        self.ax.scatter(x,y,s=1)

    def make_kde(self,x,y,h=None):
        assert h is None or h in kde_bandwidth_types

        values=np.vstack([x,y])
        if h is None:
            kde = stats.gaussian_kde(values)
        elif h in ['silverman','silverman1986']:
            kde = stats.gaussian_kde(values,'silverman')
        elif h is 'chiu1999':
            h = Chiu1999_h(values)
            kde = stats.gaussian_kde(value,h)
        else:
            raise ValueError(h)

        return kde

    
