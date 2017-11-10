import os, time

import numpy as np
import scipy.stats
import pandas as pd
import math

from bokeh.io import output_file, show, curdoc, output_notebook
from bokeh.layouts import gridplot, row, column, widgetbox, layout
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.models.glyphs import Circle, Quad, Line, VBar
from bokeh.plotting import figure, curdoc
from bokeh.models.widgets import Select, Slider
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.transform import factor_cmap

import pypospack.visualization as visualization


'''
score += 1 for err less than mean
'''


def _square_errors(total_df, err_names):
    """
    :return: (list) list of squared err names
    """
    sq_err_names = []
    for names in err_names:
        squared_names = names+'.sq'
        sq_err_names.append(squared_names)
        total_df[squared_names] = total_df[names]**2

    return sq_err_names

def get_score(df, name_list):

    name_avg_dict = {}
    name_score_dict = {}
    for names in name_list:
        name_avg_dict[names] = (df[names].mean())**2

    for names in name_list:
        name_score_dict[names] = 0
        for i in range(len(df)):
            if (df[names].iloc[i])**2 < name_avg_dict[names]:
                name_score_dict[names] += 1

    return name_score_dict




class CreateInterface(object):

    def __init__(self, data_dir='data', filename='culled_009.out'):
        vizdemo = visualization.ParetoOptimizationParamVsErrorScatter()
        vizdemo.load_data_file(fname=os.path.join(data_dir, filename))

        self.total_pandas_df = vizdemo.total_df
        self.param_names = vizdemo.param_names
        self.err_names = vizdemo.err_names
        self.sq_err_names = _square_errors(self.total_pandas_df, self.err_names)
        self.free_param_names = ['chrg_Mg', 'MgO_A', 'MgO_rho', 'OO_A', 'OO_rho', 'OO_C']
        self.qoi_names = vizdemo.qoi_names


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
        self.slider_widget = Slider(
            start=0,
            end=1,
            step=.1,
            title='test title'
        )

    def generate_plots(self, doc):
        self.bar_chart = {}

        self.bar_chart['plot_width'] = 1200
        self.bar_chart['plot_height'] = 600
        self.bar_chart['title'] = 'Score by QOI Parameter'

        self.bar_chart['object_figure'] = figure(width=self.bar_chart['plot_width'],
                                                 height=self.bar_chart['plot_height'],
                                                 title=self.bar_chart['title'],
                                                 x_range=self.qoi_names)

        score_dict = get_score(self.total_pandas_df, self.err_names)
        scores = [v for k, v in score_dict.items()]
        scores = [(s/len(self.total_pandas_df))*100 for s in scores]
        self.bar_chart['source'] = ColumnDataSource(data=dict(height=[],
                                                              x=[],
                                                              width=[],
                                                              bottom=[]))

        self.bar_chart['source'].data = {'height': scores, 'x': self.qoi_names,
                                         'width': [0.5 for i in range(len(scores))],
                                         'bottom': [0 for i in range(len(scores))]}

        self.bar_chart['bar_glyph'] = VBar(x='x', width='width', bottom='bottom', top='height', fill_color='#3264C8')

        self.bar_chart['object_figure'].xaxis.major_label_orientation = 'vertical'
        self.bar_chart['object_figure'].yaxis.axis_label = 'Percent of Simulations Below Error Threshold'

        self.bar_chart['object_figure'].add_glyph(self.bar_chart['source'], self.bar_chart['bar_glyph'])

        self.bar_chart['object_figure'].add_tools(HoverTool(tooltips=[("Score", "@height"), ("Name", "@x")]))

        doc.add_root(self.bar_chart['object_figure'])



if __name__ == "__main__":
    c_i = CreateInterface()
    c_i.create_widgets()
    c_i.start_bokeh_server()
