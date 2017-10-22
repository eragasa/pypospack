"""

This module contains classes for visualization of statistics on the data generated from the pyopspack package library

Author: Seaton Ullberg, 2017

Version Requirements:
    bokeh > 0.12.7, to avoid tornado conflicts, https://github.com/bokeh/bokeh/issues/6152
"""

import os, time

import numpy as np
import pandas as pd

import bokeh.layouts
from bokeh.io import output_file, show, curdoc, output_notebook
from bokeh.layouts import gridplot, row, column
from bokeh.models import ColumnDataSource
from bokeh.models.callbacks import CustomJS
from bokeh.models.glyphs import Circle, Quad, Line
from bokeh.plotting import figure, curdoc
from bokeh.client import push_session
from bokeh.models.widgets import Select, PreText, TextInput
from bokeh.models import Range1d
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application

import pypospack.visualization as visualization

import pylab as plb
import matplotlib.pyplot as plt
import scipy.special as special
import statistics


class Statistics_1D(object):

    def __init__(self):
        self.total_df = vizdemo.total_df
        self.param_names = vizdemo.param_names
        self.err_names = vizdemo.err_names

    def analyze_params(self, val):
        mu = statistics.mean(self.total_df[val])
        sigma = statistics.stdev(self.total_df[val])
        self.param_mu = mu
        self.param_sigma = sigma

    def analyze_err(self):
        pass


    def nix(self, val, lst):
        return [x for x in lst if x != val]

    def gaussian_fit_params(self, data, mu, sigma):
        hist, edges = np.histogram(data, density=True, bins=50)
        x = np.linspace(min(self.param_graph['plotting_data']), max(self.param_graph['plotting_data']), 1000)

        pdf = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
        #cdf = (1 + special.erf((x - mu) / np.sqrt(2 * sigma ** 2))) / 2

        self.param_hist = hist
        self.param_edges = edges
        self.param_x = x
        self.param_pdf = pdf
        #self.param_cdf = cdf

    def plot_data(self, doc):
        output_file('param_err_stats.html')
        '''
        -----------------------------------------------------------------------------------
        Define Param 
        -----------------------------------------------------------------------------------
        '''
        self.param_graph = {}
        self.param_graph['data_select'] = Select(
            title='Param Data Selection',
            value=self.param_names[0],
            options=self.param_names
        )

        self.param_graph['plotting_data'] = np.array(self.total_df[self.param_graph['data_select'].value])

        self.param_graph['plot_width'] = 610
        self.param_graph['plot_height'] = 400
        self.param_graph['title'] = 'Gaussian Fit of '+self.param_graph['data_select'].value
        self.param_graph['object_figure'] = figure(width=self.param_graph['plot_width'],
                                                   height=self.param_graph['plot_height'],
                                                   title=self.param_graph['title'])
        self.param_graph['object_figure'].min_border_left = 75
        self.param_graph_hist_source = ColumnDataSource(data=dict(
            hist=[],
            left_edge=[],
            right_edge=[]
        ))
        self.param_graph_line_source = ColumnDataSource(data=dict(
            x=[],
            y_pdf=[]
        ))

        self.analyze_params(self.param_graph['data_select'].value)
        self.gaussian_fit_params(self.total_df[self.param_graph['data_select'].value], self.param_mu, self.param_sigma)

        self.param_graph_hist_source.data = {'hist':self.param_hist, 'left_edge':self.param_edges[:-1],
                                             'right_edge':self.param_edges[1:]}
        self.param_graph_line_source.data = {'x': self.param_x,
                                             'y_pdf': self.param_pdf}


        self.param_graph['quad_glyph'] = Quad(top='hist', bottom=0, left='left_edge', right='right_edge')
        self.param_graph['pdf_glyph'] = Line(x='x', y='y_pdf', line_color="#D95B43", line_width=8, line_alpha=0.7)


        self.param_graph['object_figure'].add_glyph(self.param_graph_hist_source, self.param_graph['quad_glyph'])
        self.param_graph['object_figure'].add_glyph(self.param_graph_line_source, self.param_graph['pdf_glyph'])


        layout = column(self.param_graph['object_figure'], self.param_graph['data_select'])
        doc.add_root(layout)

        def param_select_callback(attrname, old, new):
            self.analyze_params(new)
            self.param_graph['plotting_data'] = np.array(self.total_df[new])
            self.gaussian_fit_params(np.array(self.total_df[new]), self.param_mu, self.param_sigma)
            self.param_graph_hist_source.data = {'hist':self.param_hist, 'left_edge':self.param_edges[:-1],
                                             'right_edge':self.param_edges[1:]}
            self.param_graph_line_source.data = self.param_graph_line_source.data = {'x': self.param_x,
                                             'y_pdf': self.param_pdf}
            self.param_graph['object_figure'].title.text = 'Gaussian Fit of ' + new

        self.param_graph['data_select'].on_change('value', param_select_callback)

    def start_bokeh_server(self):
        self.bokeh_app = Application(
            FunctionHandler(self.plot_data))
        self.bokeh_server = Server(
            {'/': self.bokeh_app}, num_procs=1
        )
        self.bokeh_server.start()

        # start io loop for bokeh_server
        self.bokeh_server.io_loop.add_callback(
            self.bokeh_server.show, '/'
        )
        self.bokeh_server.io_loop.start()

if __name__ == "__main__":
    data_dir = 'data'
    filename = 'culled_009.out'

    vizdemo = visualization.ParetoOptimizationParamVsErrorScatter()
    vizdemo.load_data_file(fname=os.path.join(data_dir, filename))

    stats_1D = Statistics_1D()
    stats_1D.start_bokeh_server()
    stats_1D.plot_data()
