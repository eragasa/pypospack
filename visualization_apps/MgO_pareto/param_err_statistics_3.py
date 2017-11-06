"""
This module contains classes for visualization of statistics on the data generated from the pyopspack package library

Author: Seaton Ullberg, 2017

Version Requirements:
    bokeh > 0.12.7, to avoid tornado conflicts, https://github.com/bokeh/bokeh/issues/6152
"""

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

import plotly.offline
import plotly.graph_objs

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


def resize(name, plot, total_df, axis):
    '''
    :param name: name of total_df col
    :param plot: bokeh figure object to be manipulated
    :param total_df: pd.Dataframe to draw from
    :param axis: either 'x' of 'y' 
    :return: None
    '''
    bounds = trim_data(name, total_df)
    if axis == 'x':
        plot.x_range.start = bounds['upper_bound']
        plot.x_range.end = bounds['lower_bound']
    elif axis == 'y':
        plot.y_range.start = bounds['upper_bound']
        plot.y_range.end = bounds['lower_bound']

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


def set_color(name_x, name_y, total_df, C=50):

    assert len(total_df[name_x]) == len(total_df[name_y])

    heat_frame = pd.DataFrame()

    x_upper = total_df[name_x].max()
    x_lower = total_df[name_x].min()
    y_upper = total_df[name_y].max()
    y_lower = total_df[name_y].min()

    incr_x = (x_upper - x_lower)/C
    incr_y = (y_upper - y_lower)/C

    x_vals = np.arange(x_lower, x_upper, incr_x)
    y_vals = np.arange(y_lower, y_upper, incr_y)

    plotting_x_list = []
    # make a list that repeats the x series len(y) times for plotting against y
    for i in range(len(y_vals)):
        for x in x_vals:
            plotting_x_list.append(x)

    # make a list that repeats the i_th value of the y series len(y) times for plotting against x
    plotting_y_list = []
    for i in range(len(y_vals)):
        for y in range(len(y_vals)):
            plotting_y_list.append(y_vals[i])

    necessary_stats = basic_stats([name_x, name_y], total_df)
    mean_x = necessary_stats[0][name_x][0]
    mean_y = necessary_stats[1][name_y][0]

    heat_frame['rect_x'] = plotting_x_list
    heat_frame['rect_y'] = plotting_y_list
    heat_frame['rect_width'] = incr_x*(x_upper-x_lower)
    heat_frame['rect_height'] = incr_y*(y_upper-y_lower)
    heat_frame['x_diff'] = (heat_frame['rect_x'] - mean_x) ** 2
    heat_frame['y_diff'] = (heat_frame['rect_y'] - mean_y) ** 2
    heat_frame['heat_values'] = (heat_frame['x_diff'] + heat_frame['y_diff']) / 2

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

        self.DATA_TABLE_HEAT_WIDTH = 600
        self.DATA_TABLE_HIST_WIDTH = 300
        self.DATA_TABLE_COL_WIDTH = 100

        self.data_table_heat = DataTable(
            source=ColumnDataSource(fill_data_table([self.select_widget_1.value, self.select_widget_2.value],
                                                    self.total_pandas_df)),
            columns=[
                TableColumn(field="labels", title="Statistics", width=int(self.DATA_TABLE_COL_WIDTH*0.6)),
                TableColumn(field="variable_1", title=self.select_widget_1.value, width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_2", title=self.select_widget_2.value, width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1_2", title=self.select_widget_1.value+' and '+self.select_widget_2.value,
                            width=int(self.DATA_TABLE_COL_WIDTH*1.8))
            ],
            width=self.DATA_TABLE_HEAT_WIDTH
        )
        self.data_table_hist = DataTable(
            source=ColumnDataSource(fill_data_table([self.select_widget_hist.value], self.total_pandas_df)),
            columns=[
                TableColumn(field="labels", title="Statistics", width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1", title=self.select_widget_hist.value, width=self.DATA_TABLE_COL_WIDTH),
            ],
            width=self.DATA_TABLE_HIST_WIDTH
        )

        self.x_range_slider = RangeSlider(
            start=self.total_pandas_df[self.select_widget_1.value].min(),
            end=self.total_pandas_df[self.select_widget_1.value].max(),
            value=(
                self.total_pandas_df[self.select_widget_1.value].min(),
                self.total_pandas_df[self.select_widget_1.value].max()
            ),
            step=(self.total_pandas_df[self.select_widget_1.value].max() - self.total_pandas_df[
                self.select_widget_1.value].min())/100
        )

        self.y_range_slider = RangeSlider(
            start=self.total_pandas_df[self.select_widget_2.value].min(),
            end=self.total_pandas_df[self.select_widget_2.value].max(),
            value=(
                self.total_pandas_df[self.select_widget_2.value].min(),
                self.total_pandas_df[self.select_widget_2.value].max()
            ),
            step=(self.total_pandas_df[self.select_widget_2.value].max() - self.total_pandas_df[
                self.select_widget_2.value].min()) / 100
        )

    def generate_plots(self, doc):
        '''
        :param doc: 
        :return: 
        '''

        '''
        ----------------------------------------------------------------------------------------------------------
        Define Histogram
        ----------------------------------------------------------------------------------------------------------
        '''
        self.hist_plot = {}

        self.hist_plot['plotting_data'] = np.array(self.total_pandas_df[self.select_widget_hist.value])
        self.hist_plot['plot_width'] = 610
        self.hist_plot['plot_height'] = 400
        self.hist_plot['title'] = 'Gaussian Fit of ' + self.select_widget_hist.value

        self.hist_plot['object_figure'] = figure(width=self.hist_plot['plot_width'],
                                                 height=self.hist_plot['plot_height'],
                                                 title=self.hist_plot['title'])

        self.hist_plot['source'] = ColumnDataSource(data=dict(
                hist=[],
                left_edge=[],
                right_edge=[]
            ))
        self.hist_plot['pdf_source'] = ColumnDataSource(data=dict(
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

        self.hist_plot['source'].data = {'hist': hist, 'left_edge': edges[:-1],
                                     'right_edge': edges[1:]}
        self.hist_plot['pdf_source'].data = {'x': x,
                                         'y_pdf': pdf}

        self.hist_plot['quad_glyph'] = Quad(top='hist', bottom=0, left='left_edge', right='right_edge')
        self.hist_plot['pdf_glyph'] = Line(x='x', y='y_pdf', line_color="#D95B43", line_width=8, line_alpha=0.7)

        self.hist_plot['object_figure'].add_glyph(self.hist_plot['source'], self.hist_plot['quad_glyph'])
        self.hist_plot['object_figure'].add_glyph(self.hist_plot['pdf_source'], self.hist_plot['pdf_glyph'])

        '''
        ----------------------------------------------------------------------------------------------------------
        Define Heatmap
        ----------------------------------------------------------------------------------------------------------
        '''

        self.heat_map = {}

        heat_frame = set_color(self.select_widget_1.value, self.select_widget_2.value, self.total_pandas_df)
        self.heat_map_source = ColumnDataSource(heat_frame)

        self.points_source = ColumnDataSource(data=dict(
            x=self.total_pandas_df[self.select_widget_1.value],
            y=self.total_pandas_df[self.select_widget_2.value]
        ))

        self.heat_map['plot_width'] = 610
        self.heat_map['plot_height'] = 400
        self.heat_map['title'] = "Heat Map of "+self.select_widget_1.value+' vs. '+self.select_widget_2.value
        self.heat_map['object_figure'] = figure(width=self.heat_map['plot_width'],
                                               height=self.heat_map['plot_height'],
                                               title=self.heat_map['title'])

        colors = ['#FF0000', '#F2000D', '#E6001A', '#D90026', '#CC0033', '#BF0040', '#B2004C', '#A60059', '#990066',
                      '#8C0073', '#800080', '#73008C', '#660099', '#5900A6', '#4D00B2', '#4000BF', '#3300CC', '#2600D9',
                      '#1900E6', '#0D00F2', '#0000FF']

        # mapper is a transform required by bokeh to generate heatmaps
        mapper = LinearColorMapper(palette=colors,
                                   low=self.heat_map_source.data['heat_values'].min(),
                                   high=self.heat_map_source.data['heat_values'].max())

        self.heat_map['circle_glyph'] = Circle(x='x', y='y',
                                              size=1, line_color=None,
                                              fill_color='black')
        self.heat_map['rect_glyph'] = Rect(x="rect_x", y="rect_y", width=heat_frame['rect_width'].loc[0],
                                          height=heat_frame['rect_width'].loc[0],
                                          fill_color={'field': 'heat_values', 'transform': mapper},
                                          line_color=None,
                                          fill_alpha=0.75)
        self.heat_map['object_figure'].add_glyph(self.points_source, self.heat_map['circle_glyph'])
        self.heat_map['object_figure'].add_glyph(self.heat_map_source, self.heat_map['rect_glyph'])

        self.heat_map['color_bar'] = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                                 ticker=BasicTicker(desired_num_ticks=len(colors)),
                                 label_standoff=6, border_line_color=None, location=(0, 0),
                                    background_fill_alpha=0.75)
        self.heat_map['object_figure'].add_layout(self.heat_map['color_bar'], 'left')

        WIDGETBOX_WIDTH = 610

        self.layout = layout([self.select_widget_1, self.select_widget_2, self.select_widget_hist],
                            [self.heat_map['object_figure'], self.hist_plot['object_figure']],
                             [self.x_range_slider, self.y_range_slider],
                            [widgetbox(self.data_table_heat, width=WIDGETBOX_WIDTH), widgetbox(self.data_table_hist, width=WIDGETBOX_WIDTH)])

        # unavoidable resize error forces two figures to be displayed at once rather that switch between
        # data can be updated easily but plot type should remain static
        # https://groups.google.com/a/continuum.io/forum/#!topic/bokeh/P2WPJ9an7IQ

        doc.add_root(self.layout)

        """
        --------------------------------------------------------------------
        Callbacks
        --------------------------------------------------------------------
        """

        def select_1_callback(attrname, old, new):

            self.points_source.data['x'] = self.total_pandas_df[new]

            self.data_table_heat.columns = [
                TableColumn(field="labels", title="Statistics", width=int(self.DATA_TABLE_COL_WIDTH * 0.6)),
                TableColumn(field="variable_1", title=new, width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_2", title=self.select_widget_2.value, width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1_2",
                            title=new + ' and ' + self.select_widget_2.value,
                            width=int(self.DATA_TABLE_COL_WIDTH * 1.8))
            ]

            heat_frame = set_color(new, self.select_widget_2.value, self.total_pandas_df)

            self.heat_map_source.data['rect_x'] = heat_frame['rect_x']
            self.heat_map_source.data['rect_width'] = heat_frame['rect_width']
            self.heat_map_source.data['heat_values'] = heat_frame['heat_values']

            self.data_table_heat.source = ColumnDataSource(fill_data_table([new, self.select_widget_2.value],
                                                                           self.total_pandas_df))

            self.heat_map['object_figure'].title.text = "Heat Map of " + new + ' vs. ' + self.select_widget_2.value

            self.x_range_slider.start = self.total_pandas_df[new].min()
            self.x_range_slider.end = self.total_pandas_df[new].max()
            self.x_range_slider.step = (self.total_pandas_df[new].max() - self.total_pandas_df[new].min())/100
            self.x_range_slider.value = (self.total_pandas_df[new].min(), self.total_pandas_df[new].max())

        self.select_widget_1.on_change('value', select_1_callback)

        def select_2_callback(attrname, old, new):

            self.points_source.data['y'] = self.total_pandas_df[new]

            self.data_table_heat.columns = [
                TableColumn(field="labels", title="Statistics", width=int(self.DATA_TABLE_COL_WIDTH * 0.6)),
                TableColumn(field="variable_1", title=self.select_widget_1.value, width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_2", title=new, width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1_2",
                            title=self.select_widget_1.value + ' and ' + new,
                            width=int(self.DATA_TABLE_COL_WIDTH * 1.8))
            ]

            heat_frame = set_color(self.select_widget_1.value, new, self.total_pandas_df)

            self.heat_map_source.data['rect_y'] = heat_frame['rect_y']
            self.heat_map_source.data['rect_height'] = heat_frame['rect_height']
            self.heat_map_source.data['heat_values'] = heat_frame['heat_values']

            self.data_table_heat.source = ColumnDataSource(fill_data_table([new, self.select_widget_2.value],
                                                                           self.total_pandas_df))

            self.heat_map['object_figure'].title.text = "Heat Map of " + self.select_widget_1.value + ' vs. ' + new

            self.y_range_slider.start = self.total_pandas_df[new].min()
            self.y_range_slider.end = self.total_pandas_df[new].max()
            self.y_range_slider.step = (self.total_pandas_df[new].max() - self.total_pandas_df[new].min()) / 100
            self.y_range_slider.value = (self.total_pandas_df[new].min(), self.total_pandas_df[new].max())

        self.select_widget_2.on_change('value', select_2_callback)

        def select_hist_callback(attrname, old, new):
            # get stats and gaus fit
            hist_stats_list = basic_stats([new], self.total_pandas_df)
            mu = hist_stats_list[0][new][0]
            sigma = hist_stats_list[0][new][2]
            gauss_dict = gaussian_fit(new, self.total_pandas_df)
            hist = gauss_dict['hist']
            edges = gauss_dict['edges']
            x = gauss_dict['x']
            self.hist_plot['source'].data = {'hist': hist, 'left_edge': edges[:-1],
                                             'right_edge': edges[1:]}

            if new not in self.sq_err_names:
                pdf = gauss_dict['pdf']
            else:
                pdf = chi_squared(new, self.total_pandas_df)

            self.hist_plot['pdf_source'].data = {'x': x,
                                                'y_pdf': pdf}

            self.hist_plot['object_figure'].title.text = "Gaussian Fit of " + new

            self.data_table_hist.source = ColumnDataSource(fill_data_table([new], self.total_pandas_df))
            self.data_table_hist.columns = [
                TableColumn(field="labels", title="Statistics", width=self.DATA_TABLE_COL_WIDTH),
                TableColumn(field="variable_1", title=new, width=self.DATA_TABLE_COL_WIDTH)
            ]

        self.select_widget_hist.on_change('value', select_hist_callback)

        def x_range_callback(attrname, old, new):
            self.heat_map['object_figure'].x_range.start = new[0]
            self.heat_map['object_figure'].x_range.end = new[1]


        def y_range_callback(attrname, old, new):
            self.heat_map['object_figure'].y_range.start = new[0]
            self.heat_map['object_figure'].y_range.end = new[1]

        self.x_range_slider.on_change('value', x_range_callback)
        self.y_range_slider.on_change('value', y_range_callback)

if __name__ == "__main__":

    plot_stats = PlotStatistics()
    plot_stats.create_widgets()
    plot_stats.start_bokeh_server()
