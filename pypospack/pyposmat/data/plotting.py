import matplotlib.pyplot as plt
from pypospack.pyposmat.data import BasePipeSegment


class PyposmatPlotter(BasePipeSegment):

    def __init__():
        super().__init__()

    def plot_by_cluster(self, x_axis, y_axis, filename):
        fig, ax = plt.subplots(111)
        for cid in set(self.df['cluster_id']):
            cluster_df = self.df.loc[self.df['cluster_id'] == cid]
            x_pts = cluster_df[x_axis]
            y_pts = cluster_df[y_axis]
            ax.scatter(x=x_pts, y=y_pts, s=1)
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)
        plt.savefig(filename)
