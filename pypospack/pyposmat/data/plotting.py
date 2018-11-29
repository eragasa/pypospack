import matplotlib
matplotlib.pyplot.switch_backend('agg')
import matplotlib.pyplot as plt
from pypospack.pyposmat.data import BasePipeSegment
from scipy.stats import kde
import numpy as np


class PyposmatPlotter(BasePipeSegment):

    def __init__(self):
        super().__init__()

    def plot_by_cluster(self, x_axis, y_axis, filename):
        fig, ax = plt.subplots()
        for cid in set(self.df['cluster_id']):
            cluster_df = self.df.loc[self.df['cluster_id'] == cid]
            x_pts = cluster_df[x_axis]
            y_pts = cluster_df[y_axis]
            ax.scatter(x=x_pts, y=y_pts, s=1)
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)
        plt.savefig(filename)

    # https://python-graph-gallery.com/85-density-plot-with-matplotlib/
    def plot_kde_by_cluster(self, x_axis, y_axis, filename):
        nbins = 250
        nclusters = len(set(self.df['cluster_id']))
        fig, ax = plt.subplots(nrows=1, ncols=nclusters)
        ax = {i: ax[i] for i in set(self.df['cluster_id'])}
        for cid in set(self.df['cluster_id']):
            cluster_df = self.df.loc[self.df['cluster_id'] == cid]
            x_pts = cluster_df[x_axis]
            y_pts = cluster_df[y_axis]
            k = kde.gaussian_kde([x_pts, y_pts])
            xi, yi = np.mgrid[x_pts.min():x_pts.max():nbins*1j, y_pts.min():y_pts.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            ax[cid].pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.get_cmap("cool"))
            ax[cid].set_xlabel(x_axis)
            ax[cid].set_ylabel(y_axis)
        plt.savefig(filename)
