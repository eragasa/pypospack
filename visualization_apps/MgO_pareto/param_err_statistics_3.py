"""
This module contains classes for visualization of statistics on the data generated from the pyopspack package library

Author: Seaton Ullberg, 2017

Version Requirements:
    bokeh > 0.12.7, to avoid tornado conflicts, https://github.com/bokeh/bokeh/issues/6152
"""

'''
TODO:
- nan statistics in sq err
- make color gradient more defined
- investigate better trim method
- figure out deg_f for chi
- figure out covariance position
'''

import os, time

import numpy as np
import scipy.stats
import pandas as pd
import math

from bokeh.io import output_file, show, curdoc, output_notebook
from bokeh.layouts import gridplot, row, column, widgetbox, layout
from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar, BasicTicker, PrintfTickFormatter
from bokeh.models.glyphs import Circle, Quad, Line
from bokeh.plotting import figure, curdoc
from bokeh.models.widgets import Select, DataTable, RadioGroup, RangeSlider, TableColumn
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.transform import transform

import pypospack.visualization as visualization

# frames and names loaded in class to prevent the use of global variables
# (pandas.DataFrame) total_df is used in all functions and is always the pandas.DataFrame which conatains all parameter and err data
# (list) param_names refer to the total_df headers which are used to reference param columns
# (list) err_names refer to the total_df headers which are used to refernce err columns

def _square_errors(total_df, param_names, err_names):
    """
    :return: (tuple) first value = list of all total_df headers (including squared errors)
                        second value is list of only squared err names
    """
    sq_err_names = []
    for names in err_names:
        squared_names = names+'.sq'
        sq_err_names.append(squared_names)
        total_df[squared_names] = total_df[names]**2

    all_names = param_names+err_names+sq_err_names
    return (all_names, sq_err_names)


def trim_data(name, total_df, z=2.58):
    """
    :param name: (str) total_df column name
    :param z: (float) z-statistic for 99% confidence level
    :return: (dict) contains keys 'upper_bound', 'lower_bound' which correspond to 
                upper and lower bounds of plotting axis
    """
    # name is the name of a total_df column
    sample_size = len(total_df[name])
    mean = np.mean(total_df[name])
    stdev = np.std(total_df[name])
    std_err = stdev / np.sqrt(sample_size)
    margin_err = z * std_err
    bounds = {}
    upper_bound = mean+margin_err
    lower_bound = mean-margin_err
    bounds['upper_bound'] = upper_bound
    bounds['lower_bound'] = lower_bound
    return bounds


def basic_stats(name_list, total_df):
    """
    :param name_list: (list) contains names of total_df columns to be analyzed
    :return: (list) list of dicts key=total_df column name, value=list of statistics about that col
                the order of the stats list is [mu, median, sigma, skew, kurtosis]
    """

    aggregate_stats_list = []

    for names in name_list:
        names_stats_dict = {}
        mu = np.mean(total_df[names])
        median = np.median(total_df[names])
        sigma = np.std(total_df[names])
        skew = total_df[names].skew()
        kurtosis = total_df[names].kurt()
        names_stats_dict[names] = [mu, median, sigma, skew, kurtosis]
        aggregate_stats_list.append(names_stats_dict)

    return aggregate_stats_list


def correlation_and_covariance(name1, name2, total_df):
    """
    :param name1: (str) name of total_df col to be analyzed
    :param name2: (str) name of total_df col to be analyzed
    :return: (dict) key='correlation' value=float(corr)
                    key='covariance' value=float(cov)
    """

    corr = np.corrcoef(total_df[name1], total_df[name2])
    corr = corr[0][1]
    cov = np.cov(total_df[name1], total_df[name2])
    cov = cov[0][1]

    return {'correlation':corr, 'covariance':cov}


def fill_data_table(name_list, total_df):

    assert type(name_list) == list
    assert len(name_list) <= 2

    data_table_frame = pd.DataFrame()

    if len(name_list) == 1:
        name = name_list[0]
        stats_info_list = basic_stats(name_list, total_df)
        stats_info_list = stats_info_list[0][name]
        mean = stats_info_list[0]
        median = stats_info_list[1]
        stdev = stats_info_list[2]
        skew = stats_info_list[3]
        kurtosis = stats_info_list[4]

        data_table_frame['labels'] = ['Mean', 'Median', 'StDev', 'Skew', 'Kurtosis']
        data_table_frame['variable_1'] = [mean, median, stdev, skew, kurtosis]
    else:
        data_table_frame['labels'] = ['Mean', 'Median', 'StDev', 'Skew', 'Kurtosis', 'Correlation', 'Covariance']
        stats_info_list = basic_stats(name_list, total_df)
        # tried to do this in a loop but got a scalar value error. maybe come back later
        mean, mean2 = stats_info_list[0][name_list[0]][0], stats_info_list[1][name_list[1]][0]

        median, median2 = stats_info_list[0][name_list[0]][1], stats_info_list[1][name_list[1]][1]

        stdev, stdev2 = stats_info_list[0][name_list[0]][2], stats_info_list[1][name_list[1]][2]

        skew, skew2 = stats_info_list[0][name_list[0]][3], stats_info_list[1][name_list[1]][3]

        kurtosis, kurtosis2 = stats_info_list[0][name_list[0]][4], stats_info_list[1][name_list[1]][4]

        data_table_frame['variable_1'] = [mean, median, stdev, skew, kurtosis, '', '']
        data_table_frame['variable_2'] = [mean2, median2, stdev2, skew2, kurtosis2, '', '']
        cor_cov_dict = correlation_and_covariance(name_list[0], name_list[1], total_df)
        data_table_frame['variable_1_2'] = ['', '', '', '', '', cor_cov_dict['correlation'], cor_cov_dict['covariance']]

    return data_table_frame


def set_color(name1, name2, total_df):
    """
    :param name1: (str) name of total_df col to be analyzed
    :param name2: (str) name of total_df col to be analyzed
    :return: (pandas.Series) colum of data containing the percentile of each value to be used when
                determining color
    """

    color_frame = pd.DataFrame()
    color_frame['raw_x'] = total_df[name1]
    color_frame['raw_y'] = total_df[name2]
    color_frame['mean_x'] = color_frame['raw_x'].mean()
    color_frame['mean_y'] = color_frame['raw_y'].mean()
    color_frame['diff_x'] = (color_frame['raw_x'] - color_frame['mean_x']) ** 2
    color_frame['diff_y'] = (color_frame['raw_y'] - color_frame['mean_y']) ** 2
    color_frame['combined_value'] = color_frame[['diff_x', 'diff_y']].mean(axis=1)
    color_frame['percentile'] = color_frame['combined_value'] / color_frame['combined_value'].max()

    return color_frame['percentile']


def gaussian_fit(name, total_df):
    """
    :param name: (str) name of total_df column to be fit
    :return: (dict) contains hist, edges, x, pdf
    """

    gaussian_dict = {}
    minimum = total_df[name].min()
    maximum = total_df[name].max()
    mu = basic_stats([name])
    mu = mu[0][name][0]
    sigma = basic_stats([name])
    sigma = sigma[0][name][2]
    hist, edges = np.histogram(total_df[name], density=True, bins=50, range=(minimum, maximum))
    x = np.linspace(minimum, maximum, len(total_df[name]))
    # probability density function
    pdf = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
    pdf = pdf[np.logical_not(np.isnan(pdf))]
    gaussian_dict['hist'] = hist
    gaussian_dict['edges'] = edges
    gaussian_dict['x'] = x
    gaussian_dict['pdf'] = pdf

    return gaussian_dict


def chi_squared(name, total_df, df=1):
    """
    :param name: (str) name of total_df column to be fit
    :param df: (int) degrees of freedom
    :return: (numpy.ndarray) values in chi squared curve
    """

    chi_pdf = scipy.stats.chi2.pdf(total_df[name], df)
    chi_pdf = chi_pdf[np.logical_not(np.isnan(chi_pdf))]

    return chi_pdf



class PlotStatistics(object):

    def __init__(self, data_dir='data', filename='culled_009.out'):

        vizdemo = visualization.ParetoOptimizationParamVsErrorScatter()
        vizdemo.load_data_file(fname=os.path.join(data_dir, filename))

        self.total_pandas_df = vizdemo.total_df
        self.param_names = vizdemo.param_names
        self.err_names = vizdemo.err_names
        name_tuple = _square_errors(self.total_pandas_df, self.param_names, self.err_names)
        self.all_names = name_tuple[0]
        self.sq_err_names = name_tuple[1]


    def start_bokeh_server(self):
        self.bokeh_app = Application(
            FunctionHandler(self.generate_plots))
        self.bokeh_server = Server(
            {'/': self.bokeh_app}, num_procs=1
        )
        self.bokeh_server.start()

        # start io loop for bokeh_server
        self.bokeh_server.io_loop.add_callback(
            self.bokeh_server.show, '/'
        )
        self.bokeh_server.io_loop.start()

    def create_widgets(self):
        self.select_widget_1 = Select(
            title='Variable 1 Selection',
            options=self.all_names,
            value=self.param_names[0]
        )
        self.select_widget_2 = Select(
            title='Variable 2 Selection',
            options=self.all_names,
            value=self.param_names[1]
        )
        self.select_widget_hist = Select(
            title='Histogram Variable Selection',
            options=self.all_names,
            value=self.all_names[0]
        )
        self.data_table_heat = DataTable(
            source=ColumnDataSource(),
            columns=[
                TableColumn(field="labels", title="Statistics"),
                TableColumn(field="Variable 1", title="Variable 1"),
                TableColumn(field="Variable 2", title=" Variable 2"),
                TableColumn(field="Var 1 and Var 2", title="Var 1 and Var 2")
            ]
        )
        self.data_table_hist = DataTable(
            source=ColumnDataSource(data_handler.general_stats_data),
            columns=[
                TableColumn(field="labels", title="Statistics"),
                TableColumn(field="Variable 1", title="Variable 1"),
                TableColumn(field="Variable 2", title=" Variable 2"),
                TableColumn(field="Var 1 and Var 2", title="Var 1 and Var 2")
            ]
        )

if __name__ == "__main__":

    plot_stats = PlotStatistics()
