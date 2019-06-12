import os
from collections import OrderedDict
import numpy as np
from numpy import linalg
import pandas as pd
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

import matplotlib.pyplot as plt
def plot_corr(df,plot_fn,size=10):
    '''Function plots a graphical correlation matrix for each pair of columns in the dataframe.

    Input:
        df: pandas DataFrame
        size: vertical and horizontal size of the plot'''

    corr = df.corr()
    fig, ax = plt.subplots(figsize=(size, size))
    ax.matshow(corr)
    plt.xticks(range(len(corr.columns)), corr.columns)
    plt.yticks(range(len(corr.columns)), corr.columns)
   
    plt.tight_layout()
    fig.savefig(plot_fn)
    plt.close(fig)

def covariance_analysis(config,data_directory,output_directory):
    assert isinstance(config,str) or isinstance(config,PyposmatConfigurationFile)
    assert isinstance(data_directory,str)

    print('data_directory:{}'.format(data_directory))
    print('output_directory:{}'.format(output_directory))

    if isinstance(config,str):
        o_config = PyposmatConfigurationFile()
        o_config.read(filename=config_fn)
    else:
        o_config = config

    for i in range(o_config.n_iterations):
        data_results_fn = os.path.join(data_directory,'pyposmat.results.{}.out'.format(i))
        data_kde_fn = os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1))

        data_results = PyposmatDataFile()
        data_results.read(filename=data_results_fn)

        data_kde = PyposmatDataFile()
        data_kde.read(filename=data_kde_fn)

        fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(10,10))


        corr = data_results.df[o_config.free_parameter_names].corr(method='pearson').values
        xtick_labels = [o_config.latex_labels[k]['label'] for k in o_config.free_parameter_names]
        ytick_labels = [o_config.latex_labels[k]['label'] for k in o_config.free_parameter_names]
        ax[0,0].matshow(corr,cmap='coolwarm',vmin=-1,vmax=1)
        for m in range(len(o_config.free_parameter_names)):
            for n in range(len(o_config.free_parameter_names)):
                ax[0,0].text(m,n,'{:.2f}'.format(corr[m,n]),color='w',ha='center')
        plt.sca(ax[0,0])
        plt.xticks(range(len(o_config.free_parameter_names)),xtick_labels)
        plt.yticks(range(len(o_config.free_parameter_names)),ytick_labels)

        corr = data_results.df[o_config.qoi_names].corr().values
        ax[0,1].matshow(corr,cmap='coolwarm',vmin=-1,vmax=1)
        for m in range(len(o_config.qoi_names)):
            for n in range(len(o_config.qoi_names)):
                ax[0,1].text(m,n,'{:.2f}'.format(corr[m,n]),color='w',ha='center')
        plt.sca(ax[0,1])
        xtick_labels = [o_config.latex_labels[k]['label'] for k in o_config.qoi_names]
        ytick_labels = [o_config.latex_labels[k]['label'] for k in o_config.qoi_names]
        plt.xticks(range(len(o_config.qoi_names)),xtick_labels)
        plt.yticks(range(len(o_config.qoi_names)),ytick_labels)

        corr = data_kde.df[o_config.free_parameter_names].corr().values
        xtick_labels = [o_config.latex_labels[k]['label'] for k in o_config.free_parameter_names]
        ytick_labels = [o_config.latex_labels[k]['label'] for k in o_config.free_parameter_names]
        ax[1,0].matshow(corr,cmap='coolwarm',vmin=-1,vmax=1)
        for m in range(len(o_config.free_parameter_names)):
            for n in range(len(o_config.free_parameter_names)):
                ax[1,0].text(m,n,'{:.2f}'.format(corr[m,n]),color='w',ha='center')
        plt.sca(ax[1,0])
        plt.xticks(range(len(o_config.free_parameter_names)),xtick_labels)
        plt.yticks(range(len(o_config.free_parameter_names)),ytick_labels)

        corr = data_kde.df[o_config.qoi_names].corr().values
        xtick_labels = [o_config.latex_labels[k]['label'] for k in o_config.qoi_names]
        ytick_labels = [o_config.latex_labels[k]['label'] for k in o_config.qoi_names]
        ax[1,1].matshow(corr,cmap='coolwarm',vmin=-1,vmax=1)
        for m in range(len(o_config.qoi_names)):
            for n in range(len(o_config.qoi_names)):
                ax[1,1].text(m,n,'{:.2f}'.format(corr[m,n]),color='w',ha='center')
        plt.sca(ax[1,1])
        plt.xticks(range(len(o_config.qoi_names)),xtick_labels)
        plt.yticks(range(len(o_config.qoi_names)),ytick_labels)

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        plot_fn = 'fig_cov_{}.png'.format(i)
        fig.savefig(os.path.join(output_directory,plot_fn))
        plt.close(fig)
        print(plot_fn)

if __name__ == "__main__":
    data_directory = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_unconstrained')
    output_directory = 'covariance_matrix_plot' 
    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    print("configuration file path:{}".format(config_fn))
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)
    
    covariance_analysis(config=o_config,
                        data_directory=data_directory,
                        output_directory=output_directory)
