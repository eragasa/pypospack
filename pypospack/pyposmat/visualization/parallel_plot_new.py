import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


class PyposmatParallelCoordinatesPlot(object):
    """Generalized parallel coordinates plotting object.

    Args:
        excluded_names (optional) (list): Names to exclude from addition to plot.

    Attributes:
        fig (matplotlib.Figure)
        axes (array of matplotlib.axes.Axes): Plotting axes for each dimension.
    """

    def __init__(self, excluded_names=[]):
        # public attributes
        self.fig = None
        self.axes = None
        # private attributes
        self._excluded_names = excluded_names
        self._is_exist = False  # set to True by self._init_plot()
        self._ncols = None  # number of data columns determined by first addition
        self._x = None  # will be set to the range of x
        self._patches = []  # store matplotlib patch objects here to generate legend at end 

    def add_dataframe(self, color, label, obj, names=None):
        """Add data to the plot directly from a pandas DataFrame.
        
        Args:
            color (obj): Color to assign to this data.
            - refer to https://matplotlib.org/users/colors.html for valid inputs.
            label (str): Identifier to include in the plot's legend.
            obj (pandas.DataFrame): DataFrame object to pull data from.
            names (optional) (list of str): Column names to extract from the dataframe.
            - Defaults to all (except excluded) if None.
        """
        if names:
            col_names = [name for name in names if name not in self._excluded_names]
        else:
            col_names = [name for name in list(obj) if name not in self._excluded_names]
        ncols = len(col_names)
        if not self._is_exist:
            self._init_plot(ncols)
        else:
            if ncols != self._ncols:
                raise ValueError("additional data must have {} columns of non-excluded data".format(self._ncols))
        obj = obj[col_names]  # remove undesirable columns
        self._add_data(color=color, label=label, df=obj)


    def add_datafile(self, color, label, obj, names=None):
        """Add data to the plot from a PyposmatDataFile.
        
        Args:
            color (obj): Color to assign to this data.
            - refer to https://matplotlib.org/users/colors.html for valid inputs.
            label (str): Identifier to include in the plot's legend.
            obj (PyposmatDataFile): PyposmatDataFile object to pull data from.
            names (optional) (list of str): Column names to extract from the dataframe.
            - Defaults to all (except excluded) if None.
        """
        obj = obj.df
        self.add_dataframe(color, label, obj, names)


    def add_reference_data(self, color, label, obj=None, names=None):
        raise NotImplementedError

    def make_plot(self, filename, xlabels, ylabel, title, 
                  ylim=None, plot_origin_line=True, legend_loc="upper right"):
        """Finalizes the plot construction.
        
        Args:
            filename (str): File name to save the plot as.
            xlabels (list of str): Labels on the x axis
            - It is critically important that the order of these labels
              has been maintained and matches that of each data addition.
            ylabel (str): Label for the y axis.
            title (str): Plot title.
            ylim (optional) (tuple): Range of the viewport (min, max) formatted.
            plot_origin_line (optional) (bool): Plot black line across 0 if True.
            legend_loc (optional) (str): Determines the location of the legend.
        """
        
        for i, ax in enumerate(self.axes):
            ax.set_xticks([self._x[i]], minor=False)  # modify xticks
            ax.set_xticklabels([xlabels[i]], rotation=80)  # modify xticks
            # change vertical view if desired
            if ylim:
                ax.set_ylim(ylim)
            # plot origin if desired
            if plot_origin_line:
                ax.plot(self._x, [0 for _ in self._x], color="black")
        # last tick requires separate addition
        self.axes[-1].set_xticks(self._x[-2:], minor=False)
        self.axes[-1].set_xticklabels(xlabels[-2:], rotation=80)

        # set the ylabel
        self.axes[0].set_ylabel(ylabel)
        # set the legend
        plt.legend(handles=self._patches, loc=legend_loc)
        # set the title
        self.fig.suptitle(title)
        # stack the subplots together
        plt.subplots_adjust(wspace=0, bottom=0.35)
        plt.savefig(filename)

    def _add_data(self, color, label, df):
        """Plots new data on the existing figure.
        
        Args:
            color (obj): Color of all additional data lines.
            label (str): Label to assign to the data lines.
            df (pandas.DataFrame): New data to plot.
            - assumes all filtering and size checking has been done.
        """
        # store the label
        patch = mpatches.Patch(color=color, label=label)
        self._patches.append(patch)

        # iterate over axes
        for i, ax in enumerate(self.axes):
            # set static x limit
            ax.set_xlim((self._x[i], self._x[i+1]))
            # iterate through the data
            for i_row, row in df.iterrows():
                ax.plot(self._x, row, color=color)

    def _init_plot(self, ncols):
        """Initialize a matplotlib plotting object.
        
        Args:
            ncols (int): Number of columns in the dataset.
        """
        self._ncols = ncols
        self._x = range(ncols)
        self.fig, self.axes = plt.subplots(nrows=1, ncols=ncols-1, figsize=(8, 4), sharey=True)
        self._is_exist = True
