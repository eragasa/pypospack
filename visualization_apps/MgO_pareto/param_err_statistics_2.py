"""
This module contains classes for visualization of statistics on the data generated from the pyopspack package library

Author: Seaton Ullberg, 2017

Version Requirements:
    bokeh > 0.12.7, to avoid tornado conflicts, https://github.com/bokeh/bokeh/issues/6152
"""

'''
TODO:
Chi Square 
remove error outside 99%
--trim data function
'''

import os, time

import numpy as np
import pandas as pd

from bokeh.io import output_file, show, curdoc, output_notebook
from bokeh.layouts import gridplot, row, column, widgetbox, layout
from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar
from bokeh.models.glyphs import Circle, Quad, Line
from bokeh.plotting import figure, curdoc
from bokeh.models.widgets import Select, DataTable, RadioGroup, RangeSlider, TableColumn
from bokeh.models import Range1d
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.transform import transform

import pypospack.visualization as visualization

import scipy.special as special
import scipy.stats
import statistics

class Plot_Statistics(object):

    def __init__(self):
        pass

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
            options=data_handler.total_names,
            value=data_handler.param_names[0]
        )

        self.select_widget_2 = Select(
            title='Variable 2 Selection',
            options=data_handler.total_names,
            value=data_handler.param_names[1]
        )

        self.select_widget_hist = Select(
            title='Histogram Variable Selection',
            options=data_handler.total_names,
            value=data_handler.total_names[0]
        )

        data_handler.general_stats(self.select_widget_1.value, self.select_widget_2.value)

        self.data_table = DataTable(
            source=ColumnDataSource(data_handler.general_stats_data),
            columns=[
                TableColumn(field="labels", title="Statistics"),
                TableColumn(field="Variable 1", title="Variable 1"),
                TableColumn(field="Variable 2", title=" Variable 2"),
                TableColumn(field="Var 1 and Var 2", title="Var 1 and Var 2")
            ]
        )

    def generate_plots(self, doc):

        hist_plot = {}
        heat_map = {}

        def generate_histogram():
            nonlocal  hist_plot
            hist_plot['plotting_data'] = np.array(data_handler.total_df[self.select_widget_hist.value])
            hist_plot['plot_width'] = 610
            hist_plot['plot_height'] = 400
            hist_plot['title'] = 'Gaussian Fit of ' + self.select_widget_hist.value

            hist_plot['object_figure'] = figure(width=hist_plot['plot_width'],
                                                     height=hist_plot['plot_height'],
                                                     title=hist_plot['title'])
            #hist_plot['object_figure'].min_border_left = 75

            hist_plot['source'] = ColumnDataSource(data=dict(
                hist=[],
                left_edge=[],
                right_edge=[]
            ))
            hist_plot['pdf_source'] = ColumnDataSource(data=dict(
                x=[],
                y_pdf=[]
            ))
            data_handler.general_stats(self.select_widget_hist.value, None)
            data_handler.gaussian_fit(data_handler.total_df[self.select_widget_hist.value],
                                      data_handler.mu,
                                      data_handler.sigma,
                                      min(data_handler.total_df[self.select_widget_hist.value]),
                                      max(data_handler.total_df[self.select_widget_hist.value]))

            hist_plot['source'].data = {'hist': data_handler.hist, 'left_edge': data_handler.edges[:-1],
                                     'right_edge': data_handler.edges[1:]}
            hist_plot['pdf_source'].data = {'x': data_handler.x,
                                         'y_pdf': data_handler.pdf}
            hist_plot['quad_glyph'] = Quad(top='hist', bottom=0, left='left_edge', right='right_edge')
            hist_plot['pdf_glyph'] = Line(x='x', y='y_pdf', line_color="#D95B43", line_width=8, line_alpha=0.7)

            hist_plot['object_figure'].add_glyph(hist_plot['source'], hist_plot['quad_glyph'])
            hist_plot['object_figure'].add_glyph(hist_plot['pdf_source'], hist_plot['pdf_glyph'])

        generate_histogram()
        self.hist_plot = hist_plot

        def generate_heatmap():
            nonlocal  heat_map

            data_handler.general_stats(self.select_widget_1.value, self.select_widget_2.value)

            heat_map['plotting_data'] = data_handler.total_df
            heat_map['source'] = ColumnDataSource(data=dict(
                x=[],
                y=[],
                heat_values=[]
            ))
            heat_map['source'].data = {'x':data_handler.total_df[self.select_widget_1.value],
                               'y':data_handler.total_df[self.select_widget_2.value],
                               'heat_values':data_handler.total_df['heat_values']}

            heat_map['plot_width'] = 610
            heat_map['plot_height'] = 400
            heat_map['title'] = "Heat Map of "+self.select_widget_1.value+' vs. '+self.select_widget_2.value
            heat_map['object_figure'] = figure(width=heat_map['plot_width'],
                                               height=heat_map['plot_height'],
                                               title=heat_map['title'])
            #heat_map['object_figure'].min_border_left = 500
            colors = ["#E50008", "#CE071A", "#B70E2C", "#A0153E", "#891C51", "#722463", "#5B2B75", "#443288",
                      "#2D399A", "#1640AC", "#0048BF"]
            mapper = LinearColorMapper(palette=colors,
                                       low=heat_map['source'].data['heat_values'].min(),
                                       high=heat_map['source'].data['heat_values'].max())
            heat_map['circle_glyph'] = Circle(x='x', y='y',
                                              size=1, line_color=None, fill_color=transform('heat_values', mapper))
            heat_map['object_figure'].add_glyph(heat_map['source'], heat_map['circle_glyph'])


        generate_heatmap()
        self.heat_map = heat_map

        self.layout = layout([self.select_widget_1, self.select_widget_2, self.select_widget_hist],
                            [self.heat_map['object_figure'], self.hist_plot['object_figure']],
                            [self.data_table])

        # unavoidable resize error forces two figures to be displayed at once rather that switch between
        # https://groups.google.com/a/continuum.io/forum/#!topic/bokeh/P2WPJ9an7IQ

        doc.add_root(self.layout)

        def select_1_callback(attrname, old, new):
            data_handler.general_stats(new, self.select_widget_2.value)
            self.heat_map['source'].data['x'] = data_handler.total_df[new]
            self.heat_map['object_figure'].title.text = "Heat Map of "+new+' vs. '+self.select_widget_2.value
            #self.heat_map['object_figure'].x_range.start = data_handler.total_df[new].min()
            #self.heat_map['object_figure'].x_range.end = data_handler.total_df[new].max()
            self.data_table.source = ColumnDataSource(data_handler.general_stats_data)
            #self.layout.children[1].children[0] = self.heat_map['object_figure']

        self.select_widget_1.on_change('value', select_1_callback)

        def select_2_callback(attrname, old, new):
            data_handler.general_stats(self.select_widget_1.value, new)
            self.heat_map['source'].data['y'] = data_handler.total_df[new]
            self.heat_map['object_figure'].title.text = "Heat Map of "+self.select_widget_1.value+' vs. '+new
            #self.heat_map['object_figure'].y_range.start = data_handler.total_df[new].min()
            #self.heat_map['object_figure'].y_range.end = data_handler.total_df[new].max()
            self.data_table.source = ColumnDataSource(data_handler.general_stats_data)
            #self.layout.children[1].children[0] = self.heat_map['object_figure']

        self.select_widget_2.on_change('value', select_2_callback)

        def select_hist_callback(attrname, old, new):
            data_handler.general_stats(new, None)
            data_handler.gaussian_fit(data_handler.total_df[new],
                                    data_handler.mu,
                                    data_handler.sigma,
                                    data_handler.total_df[new].min(),
                                    data_handler.total_df[new].max())
            hist_plot['source'].data = {'hist': data_handler.hist, 'left_edge': data_handler.edges[:-1],
                                        'right_edge': data_handler.edges[1:]}

            if new in data_handler.squared_err_names:
                data_handler.chi_squared(data_handler.total_df[new])
                hist_plot['pdf_source'].data = {'x': data_handler.x,
                                            'y_pdf': data_handler.chi_pdf}
            else:
                hist_plot['pdf_source'].data = {'x': data_handler.x,
                                                'y_pdf': data_handler.pdf}
            self.hist_plot['object_figure'].title.text = 'Gaussian Fit of ' + new
            #self.layout.children[1].children[0] = self.hist_plot['object_figure']

        self.select_widget_hist.on_change('value', select_hist_callback)

class Statistics_Data_Handler(object):

    def __init__(self):
        self.total_df = vizdemo.total_df
        self.param_names = vizdemo.param_names
        self.err_names = vizdemo.err_names
        # create new columns for squared error
        self.squared_err_names = []
        for names in self.err_names:
            squared_err_name = names+'.sq'
            self.squared_err_names.append(squared_err_name)
            self.total_df[squared_err_name] = self.total_df[names]**2

        self.total_names = self.param_names + self.err_names + self.squared_err_names

        # trim values outside 99% confidence
        self.err_margin_dict = {}

        for err_names in self.err_names:
            self.trim_data(err_names)

        for sq_err_names in self.squared_err_names:
            self.trim_data(sq_err_names)

    def trim_data(self, name):
        # name is the name of a total_df column
        sample_size = len(self.total_df[name])
        mean = statistics.mean(self.total_df[name])
        try:
            stdev = statistics.stdev(self.total_df[name])
            z = 2.58
            std_err = stdev / np.sqrt(sample_size)
            margin_err = z * std_err
            self.err_margin_dict[name] = margin_err
        except OverflowError:
            self.err_margin_dict[name] = statistics.mean(self.total_df[name])

        self.total_df[name] = self.total_df[self.total_df[name].between(mean-self.err_margin_dict[name], mean+self.err_margin_dict[name])]
        self.total_df[name].dropna()

    def general_stats(self, val, val2):
        # val is the name associated with a total_df column
        # val for var1 and val2 for var2
        mu = statistics.mean(self.total_df[val])
        median = statistics.median(self.total_df[val])
        sigma = statistics.stdev(self.total_df[val])
        skew = scipy.stats.skew(self.total_df[val])
        kurtosis = scipy.stats.kurtosis(self.total_df[val])

        # copy for use
        self.mu, self.mu2 = mu, None
        self.median, self.median2 = median, None
        self.sigma, self.sigma2 = sigma, None
        self.skew, self.skew2 = skew, None
        self.kurtosis, self.kurtosis2 = kurtosis, None
        self.corr = None
        self.cov = None

        if val2 != None:
            mu2 = statistics.mean(self.total_df[val2])
            median2 = statistics.median(self.total_df[val2])
            sigma2 = statistics.stdev(self.total_df[val2])
            skew2 = scipy.stats.skew(self.total_df[val2])
            kurtosis2 = scipy.stats.kurtosis(self.total_df[val2])

            self.mu2 = mu2
            self.median2 = median2
            self.sigma2 = sigma2
            self.skew2 = skew2
            self.kurtosis2 = kurtosis2

            correlation = np.corrcoef(x=self.total_df[val], y=self.total_df[val2])
            correlation = correlation[0][1]
            covariance = np.cov(m=self.total_df[val], y=self.total_df[val2])
            covariance = covariance[0][0]
            '''
            Which covariance position do I use?
            '''
            self.corr = correlation
            self.cov = covariance
            # obtain heat map 'values'
            x_heats = []
            for values in self.total_df[val]:
                heat_x = 1/((values-self.mu)**2)
                x_heats.append(heat_x)
            y_heats = []
            for values in self.total_df[val2]:
                heat_y = 1/((values-self.mu2)**2)
                y_heats.append(heat_y)
            heat_zip = zip(x_heats, y_heats)
            avg_heat = [sum(e) / len(e) for e in heat_zip]
            self.total_df['heat_values'] = avg_heat


        self.general_stats_data = pd.DataFrame()
        self.general_stats_data['labels'] = ['Mean', 'Median', 'StDev', 'Skew', 'Kurtosis', 'Correlation', 'Covariance']
        self.general_stats_data['Variable 1'] = [self.mu, self.median, self.sigma, self.skew, self.kurtosis, '', '']
        self.general_stats_data['Variable 2'] = [self.mu2, self.median2, self.sigma2, self.skew2, self.kurtosis2, '', '']
        self.general_stats_data['Var 1 and Var 2'] = ['', '', '', '', '', self.corr, self.cov]

    def gaussian_fit(self, data, mu, sigma, minimum, maximum):
        # data is a total_df data column
        # mu is the mean of data
        # sigma is the StDev of data
        hist, edges = np.histogram(data, density=True, bins=50)
        x = np.linspace(minimum, maximum, 1000)
        # probability density function
        pdf = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
        self.hist = hist
        self.edges = edges
        self.x = x
        self.pdf = pdf

    def chi_squared(self, data):
        # data is a total_df data column
        # how many degrees of freedom do I use?
        chi_pdf = scipy.stats.chi2.pdf(data, 1)
        self.chi_pdf = chi_pdf

if __name__ == "__main__":
    data_dir = 'data'
    filename = 'culled_009.out'

    vizdemo = visualization.ParetoOptimizationParamVsErrorScatter()
    vizdemo.load_data_file(fname=os.path.join(data_dir, filename))

    data_handler = Statistics_Data_Handler()
    plot_stats = Plot_Statistics()
    plot_stats.create_widgets()
    plot_stats.start_bokeh_server()


