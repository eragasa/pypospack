import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatAbstractPlot

class PyposmatParallelCoordinatesPlot(PyposmatAbstractPlot):

    normalize_types = ['by_qoi_target']
    def __init__(self,config=None,data=None,excluded_names=[]):
        assert isinstance(excluded_names,list)

        PyposmatAbstractPlot.__init__(self,config=config,data=data)
        self.legend_patches = []
        self.x_limits = None
        self.y_limits = None
        self.excluded_names = excluded_names

    def create_subplots(self,
                        figsize=(10,5)):
        PyposmatAbstractPlot.create_subplots(self,
                nrows=1,ncols=1,
                sharex=False,sharey=False,
                figsize = figsize,
                squeeze=True,
                subplot_kw=None,gridspec_kw=None)

    def plot(self,
             config=None,data=None,
             normalize_type='by_qoi_target',
             label=None,color=None,nsmallest=20,
             linewidth=1,
             alpha=.7,
             cluster_id=None):
        """
        Args:
            config(None)(str)(PyposmatConfigurationFile)
            data(None)(str)(PyposmatDataFile)
            normalize_type(str)
            label(None)(str)
            color(None)(str)
            nsmallest=20
        """
        if self.fig is None or self.ax is None:
            self.create_subplots()

        self.initialize_configuration(config=config)
        self.initialize_data(data=data)
        self.data.create_normalized_errors(
                normalize_type='by_qoi_target',
                qoi_targets=self.configuration.qoi_targets)

        normalized_error_names = self.configuration.normalized_error_names
        self.data.df['score'] = self.data.df[normalized_error_names].abs().sum(axis=1)

        # data to plot
        x = range(len(normalized_error_names))
        if cluster_id is None:
            df = self.data.df.nsmallest(20,'score')[normalized_error_names]
        else:
            df = self.data.df.loc[self.data.df['cluster_id'] == cluster_id]
            df = df.nsmallest(20,'score')[normalized_error_names]

        if isinstance(color,str) and isinstance(label,str):
            self.legend_patches.append(mpatches.Patch(color=color,label=label))
        
        # plot the data
        for i, row in df.iterrows():
            self.ax.plot(x,
                    row,
                    color=color,
                    linewidth=linewidth,
                    alpha=alpha)

        if self.x_limits is None:
            self.set_xlimits(
                    x_lim_min=min(x),
                    x_lim_max=max(x))
        else:
            x_lim_min = min(min(x),self.x_limits[0])
            x_lim_max = max(max(x),self.x_limits[1])
            self.set_xlimits(
                    x_lim_min=x_lim_min,
                    x_lim_max=x_lim_max)

        #if self.y_limits is None:
        #    y_lim_min = -max(df.abs().max())
        #    y_lim_max = max(df.abs().max())
        #    self.set_ylimits(
        #            y_lim_min=y_lim_min,
        #            y_lim_max=y_lim_max)
        #else:
        #    y_lim_min = -max(self.y_limits[1],df.abs().max())
        #    y_lim_max = max(self.y_limits[1],df.abs().max())
        #    self.set_ylimits(
        #            y_lim_min=y_lim_min,
        #            y_lim_max=y_lim_max)


    def plot_reference_potentials(self,
            config=None,data=None,
            normalize_type='by_qoi_target',
            labels=None,colors=None,nsmallest=20,
            linewidth=2,
            alpha=1.):
        """

        Args:
            config(str)(PyposmatConfigurationFile)
            data(str)(PyposmatDataFile)
            normalize_type(str)
            labels(None)(list of str)(dict)
            colors(None((list of str)(dict)
            nsmallest(int)
            alpha(float)
        """
        assert isinstance(config,str) or isinstance(config,PyposmatConfigurationFile)
        assert isinstance(data,str) or isinstance(data,PyposmatDataFile)
        assert isinstance(normalize_type,str)
        assert labels is None \
                or isinstance(labels,list) \
                or isinstance(labels,dict)
        assert colors is None \
                or isinstance(colors,list) \
                or isinstance(colors,dict)
        assert isinstance(nsmallest,int)
        assert isinstance(alpha,float)

        self.initialize_configuration(config=config)
        self.initialize_data(data=data)

        self.data.create_normalized_errors(
                normalize_type='by_qoi_target',
                qoi_targets=self.configuration.qoi_targets)

        if normalize_type == 'by_qoi_target':
            normalized_error_names = self.configuration.normalized_error_names
            self.data.df['score'] = self.data.df[normalized_error_names].abs().sum(axis=1)

            names = self.configuration.normalized_error_names
        else:
            raise ValueError('unknown normalize_type:{}'.format(normalize_type))

        x = range(len(normalized_error_names))

        self.reference_handles = []
        for i,row in self.data.df.iterrows():
            if labels is None:
                label = row['sim_id']
            else:
                label = labels[i]

            if colors is None:
                color = None
            else:
                color = colors[i]

            self.reference_handles += self.ax.plot(
                    x,
                    row[names],
                    color=color,
                    label=label,
                    linewidth=linewidth,
                    alpha=alpha)

    def set_legend(self,handles=None,location=None):
        if handles is None:
            print(type(self.legend_patches))
            print(self.legend_patches)
            print(type(self.reference_handles))
            print(self.reference_handles)
            handles_ = self.legend_patches + self.reference_handles

        if location is None:
            location_ = 'upper right'

        plt.legend(handles=handles_,loc=location)

    def set_xlabel(self,label,ax=None):
        if ax is None:
            ax_ = self.ax
        ax_.set_xlabel(label)
    
    def set_ylabel(self,label,ax=None):
        if ax is None:
            ax_ = self.ax
        ax_.set_ylabel(label)

    def set_xticks(self,config,names,ax=None):
        if isinstance(config,PyposmatConfigurationFile):
            config_ = config
        elif isinstance(config,str):
            config_ = PyposmatConfigurationFile()
            config_.read(config)
        else:
            raise TypeError()

        if names == 'qoi_names':
            names_ = config_.qoi_names
        elif isinstance(names,list):
            names_ = names
        else:
            raise TypeError()

        if ax is None:
            ax_ = self.ax
        latex_labels = [config_.latex_labels[k]['label'] for k in names_]
        print(latex_labels)
        ax_.set_xticks(range(len(latex_labels)))
        ax_.set_xticklabels(latex_labels)

    def set_xlimits(self,x_lim_min,x_lim_max,ax=None):
        if ax is None:
            ax = self.ax
        self.x_limits = [x_lim_min,x_lim_max]
        ax.set_xlim(x_lim_min,x_lim_max)

    def set_ylimits(self,y_lim_min,y_lim_max,ax=None):
        if ax is None:
            ax = self.ax
        self.y_limits = [y_lim_min,y_lim_max]
        ax.set_ylim(y_lim_min,y_lim_max)

    def show_figure(self):
        plt.tight_layout()
        plt.show()

    def save_figure(self,filename,dpi=1300):
        assert isinstance(filename,str)
        assert isinstance(dpi,int)
        plt.tight_layout()
        self.fig.savefig(filename,dpi=dpi)

class PyposmatMultipleParallelCoordinatesPlot(PyposmatParallelCoordinatesPlot):

    def __init__(self,excluded_names=[]):
        PyposmatAbstractPlot.__init__(self)

    def plot(self,configs,datas):
        assert isinstance(configs,OrderedDict)
        assert isinstance(datas,OrderedDict)
        for k,v in configs.items():
            assert isinstance(v,PyposmatConfigurationFile) or isinstance(v,str)
        for k in datas:
            assert isinstance(v,PyposmatDataFile) or isinstance(v,str)
        assert list(configs.keys()) == list(configs.keys())


        for i in range(len(configs)):
            PyposmatParallelCoordinatesPlot(configs[i],data[i],color[i])

if __name__ == "__main__":
    import os
    import pypospack.utils
    o = PyposmatParallelCoordinatesPlot()
    o.plot(config=os.path.join(
                pypospack.utils.get_pypospack_root_directory(),
                'data','Si__sw__data','pareto_optimization_p_3.5_q_0.5',
                'pyposmat.config.in'),
           data=os.path.join(
                pypospack.utils.get_pypospack_root_directory(),
                'data','Si__sw__data','pareto_optimization_p_3.5_q_0.5',
                'pyposmat.kde.19.out'),
           color='blue')
    o.show_figure()

    o = PyposmatMultipleParallelCoordinatesPlot()
