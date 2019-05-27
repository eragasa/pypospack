import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from collections import OrderedDict

def make_pareto_plot_3d(df,idx):
    fig_height = 1
    fig_width = 2
    fig_aspect = fig_height/fig_width
    fig = plt.figure(figsize=plt.figaspect(fig_aspect))

    col_names = list(df.columns.values)

    ax=OrderedDict()
    ax['all'] = fig.add_subplot(1,2,1,projection='3d')
    ax['all'].scatter(
        df[col_names[0]],
        df[col_names[1]],
        df[col_names[2]],
        s=1)
    ax['pareto'] = fig.add_subplot(1,2,2,projection='3d')
    ax['pareto'].scatter(
        df.loc[idx,col_names[0]],
        df.loc[idx,col_names[1]],
        df.loc[idx,col_names[2]],
        s=1)
    plt.show()

def make_pareto_plot_2d(df,idx):
    fig,ax = plt.subplots(nrows=1,ncols=2)
    ax[0].scatter(
        df['x'],
        df['y'],
        s=1)
    ax[1].scatter(
        df.loc[idx,'x'],
        df.loc[idx,'y'],
        s=1)
    plt.show()
