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
from bokeh.models.glyphs import Circle, Quad, Line, Rect
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

'''
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
'''

def set_color(name_x, name_y, total_df, C=500):

    assert len(total_df[name_x]) == len(total_df[name_y])

    heat_frame = pd.DataFrame()
    '''
    x_upper = total_df[name_x].max()
    x_lower = total_df[name_x].min()
    y_upper = total_df[name_y].max()
    y_lower = total_df[name_y].min()
    '''

    x_bounds = trim_data(name_x, total_df)
    x_upper = x_bounds['upper_bound']
    x_lower = x_bounds['lower_bound']
    #x_upper = 2*x_upper
    #x_lower = 0.5*x_lower

    y_bounds = trim_data(name_y, total_df)
    y_upper = y_bounds['upper_bound']
    y_lower = y_bounds['lower_bound']
    #y_upper = 2*y_upper
    #y_lower = 0.5*y_lower

    incr_x = (x_upper - x_lower)/C
    incr_y = (y_upper - y_lower)/C

    x_vals = np.arange(x_lower, x_upper, incr_x)
    y_vals = np.arange(y_lower, y_upper, incr_y)

    necessary_stats = basic_stats([name_x, name_y], total_df)
    mean_x = necessary_stats[0][name_x][0]
    mean_y = necessary_stats[1][name_y][0]


    all_x_list = [list(x_vals) for i in range(len(y_vals))]
    all_y_list = [list(y_vals) for i in range(len(x_vals))]

    plotting_x_list = []
    for sublists_x in all_x_list:
        for elements_x in sublists_x:
            plotting_x_list.append(elements_x)

    plotting_y_list = []
    list_count = 0
    for sublists_y in all_y_list:
        if list_count != 0:
            sublists_y = sublists_y[list_count:]+sublists_y[:list_count]
        list_count += 1
        for elements_y in sublists_y:
            plotting_y_list.append(elements_y)

    heat_frame['heat_x'] = plotting_x_list
    heat_frame['heat_y'] = plotting_y_list
    heat_frame['rect_width'] = incr_x/2
    heat_frame['rect_height'] = incr_y/2
    heat_frame['x_diff'] = (heat_frame['heat_x'] - mean_x) ** 2
    heat_frame['y_diff'] = (heat_frame['heat_y'] - mean_y) ** 2
    heat_frame['avg_diff'] = (heat_frame['x_diff'] + heat_frame['y_diff']) / 2

    '''
    color_frame['x'] = x_vals
    color_frame['y'] = y_vals
    
    color_frame['x_diff'] = (color_frame['x'] - mean_x)**2
    color_frame['y_diff'] = (color_frame['y'] - mean_y)**2
    color_frame['avg_diff'] = (color_frame['x_diff'] + color_frame['y_diff'])/2
    '''

    return heat_frame


def gaussian_fit(name, total_df):
    """
    :param name: (str) name of total_df column to be fit
    :return: (dict) contains hist, edges, x, pdf
    """

    gaussian_dict = {}
    minimum = total_df[name].min()
    maximum = total_df[name].max()
    mu = basic_stats([name], total_df)
    mu = mu[0][name][0]
    sigma = basic_stats([name], total_df)
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

        DATA_TABLE_HEAT_WIDTH = 500
        DATA_TABLE_HIST_WIDTH = 300
        DATA_TABLE_COL_WIDTH = 75

        self.data_table_heat = DataTable(
            source=ColumnDataSource(fill_data_table([self.select_widget_1.value, self.select_widget_2.value],
                                                    self.total_pandas_df)),
            columns=[
                TableColumn(field="labels", title="Statistics", width=DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1", title=self.select_widget_1.value, width=DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_2", title=self.select_widget_2.value, width=DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1_2", title=self.select_widget_1.value+' and '+self.select_widget_2.value,
                            width=DATA_TABLE_COL_WIDTH)
            ],
            width=DATA_TABLE_HEAT_WIDTH
        )
        self.data_table_hist = DataTable(
            source=ColumnDataSource(fill_data_table([self.select_widget_hist.value], self.total_pandas_df)),
            columns=[
                TableColumn(field="labels", title="Statistics", width=DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1", title=self.select_widget_hist.value, width=DATA_TABLE_COL_WIDTH),
            ],
            width=DATA_TABLE_HIST_WIDTH
        )

    def generate_plots(self, doc):

        # use nested functions to get around the (doc) parameter restriction

        hist_plot = {}
        heat_map = {}

        def generate_histogram():
            nonlocal  hist_plot
            hist_plot['plotting_data'] = np.array(self.total_pandas_df[self.select_widget_hist.value])
            hist_plot['plot_width'] = 610
            hist_plot['plot_height'] = 400
            hist_plot['title'] = 'Gaussian Fit of ' + self.select_widget_hist.value

            hist_plot['object_figure'] = figure(width=hist_plot['plot_width'],
                                                     height=hist_plot['plot_height'],
                                                     title=hist_plot['title'])

            hist_plot['source'] = ColumnDataSource(data=dict(
                hist=[],
                left_edge=[],
                right_edge=[]
            ))
            hist_plot['pdf_source'] = ColumnDataSource(data=dict(
                x=[],
                y_pdf=[]
            ))

            # get stats and gaus fit
            hist_stats_list = basic_stats([self.select_widget_hist.value], self.total_pandas_df)
            mu = hist_stats_list[0][self.select_widget_hist.value][0]
            sigma = hist_stats_list[0][self.select_widget_hist.value][2]
            gauss_dict = gaussian_fit(self.select_widget_hist.value, self.total_pandas_df)
            hist = gauss_dict['hist']
            edges = gauss_dict['edges']
            x = gauss_dict['x']
            pdf = gauss_dict['pdf']

            hist_plot['source'].data = {'hist': hist, 'left_edge': edges[:-1],
                                     'right_edge': edges[1:]}
            hist_plot['pdf_source'].data = {'x': x,
                                         'y_pdf': pdf}

            hist_plot['quad_glyph'] = Quad(top='hist', bottom=0, left='left_edge', right='right_edge')
            hist_plot['pdf_glyph'] = Line(x='x', y='y_pdf', line_color="#D95B43", line_width=8, line_alpha=0.7)

            hist_plot['object_figure'].add_glyph(hist_plot['source'], hist_plot['quad_glyph'])
            hist_plot['object_figure'].add_glyph(hist_plot['pdf_source'], hist_plot['pdf_glyph'])

        generate_histogram()
        self.hist_plot = hist_plot

        def generate_heatmap():
            nonlocal  heat_map
            # get stats
            #heat_stats_list = basic_stats([self.select_widget_1.value, self.select_widget_2.value], self.total_pandas_df)

            heat_map['rect_source'] = ColumnDataSource(data=dict(
                heat_x=[],
                heat_y=[],
                heat_values=[]
            ))
            heat_map['point_source'] = ColumnDataSource(data=dict(
                x=[],
                y=[]
            ))
            # determine percentile of each value for color assignment
            #frames_list = set_color(self.select_widget_1.value, self.select_widget_2.value, self.total_pandas_df)
            heat_frame = set_color(self.select_widget_1.value, self.select_widget_2.value, self.total_pandas_df)
            heat_map['rect_source'].data = {'heat_x': heat_frame['heat_x'],
                                       'heat_y': heat_frame['heat_y'],
                                       'heat_values': heat_frame['avg_diff']}
            heat_map['point_source'].data = {'x': self.total_pandas_df[self.select_widget_1.value],
                                             'y': self.total_pandas_df[self.select_widget_2.value]
            }
            heat_map['plot_width'] = 610
            heat_map['plot_height'] = 400
            heat_map['title'] = "Heat Map of "+self.select_widget_1.value+' vs. '+self.select_widget_2.value
            heat_map['object_figure'] = figure(width=heat_map['plot_width'],
                                               height=heat_map['plot_height'],
                                               title=heat_map['title'])
            x_bounds = trim_data(self.select_widget_1.value, self.total_pandas_df)
            y_bounds = trim_data(self.select_widget_2.value, self.total_pandas_df)
            x_upper = x_bounds['upper_bound']
            x_lower = x_bounds['lower_bound']
            y_upper = y_bounds['upper_bound']
            y_lower = y_bounds['lower_bound']

            heat_map['object_figure'].x_range.start = x_lower
            heat_map['object_figure'].x_range.end = x_upper
            heat_map['object_figure'].y_range.start = y_lower
            heat_map['object_figure'].y_range.end = y_upper

            colors = ['#FF0000', '#F2000D', '#E6001A', '#D90026', '#CC0033', '#BF0040', '#B2004C', '#A60059', '#990066',
                      '#8C0073', '#800080', '#73008C', '#660099', '#5900A6', '#4D00B2', '#4000BF', '#3300CC', '#2600D9',
                      '#1900E6', '#0D00F2', '#0000FF']
            # mapper is a transform required by bokeh to generate heatmaps
            mapper = LinearColorMapper(palette=colors,
                                       low=heat_map['rect_source'].data['heat_values'].min(),
                                       high=heat_map['rect_source'].data['heat_values'].max())

            heat_map['circle_glyph'] = Circle(x='x', y='y',
                                              size=1, line_color=None,
                                              fill_color='black')
            heat_map['rect_glyph'] = Rect(x="heat_x", y="heat_y", width=heat_frame['rect_width'].loc[0],
                                          height=heat_frame['rect_width'].loc[0],
                                          fill_color={'field': 'heat_values', 'transform': mapper},
                                          line_color=None)
            heat_map['object_figure'].add_glyph(heat_map['point_source'], heat_map['circle_glyph'])
            heat_map['object_figure'].add_glyph(heat_map['rect_source'], heat_map['rect_glyph'])

            heat_map['color_bar'] = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                                 ticker=BasicTicker(desired_num_ticks=len(colors)),
                                 label_standoff=6, border_line_color=None, location=(0, 0))
            heat_map['object_figure'].add_layout(heat_map['color_bar'], 'left')

        generate_heatmap()
        self.heat_map = heat_map

        WIDGETBOX_WIDTH = 610

        self.layout = layout([self.select_widget_1, self.select_widget_2, self.select_widget_hist],
                            [self.heat_map['object_figure'], self.hist_plot['object_figure']],
                            [widgetbox(self.data_table_heat, width=WIDGETBOX_WIDTH), widgetbox(self.data_table_hist, width=WIDGETBOX_WIDTH)])

        # unavoidable resize error forces two figures to be displayed at once rather that switch between
        # data can be updated easily but plot type should remain static
        # https://groups.google.com/a/continuum.io/forum/#!topic/bokeh/P2WPJ9an7IQ

        doc.add_root(self.layout)

if __name__ == "__main__":

    plot_stats = PlotStatistics()
    plot_stats.create_widgets()
    plot_stats.start_bokeh_server()
    # print(set_color(plot_stats.param_names[0], plot_stats.param_names[1], plot_stats.total_pandas_df))
