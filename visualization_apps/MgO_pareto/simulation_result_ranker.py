import os

import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.models.glyphs import Quad, Line
from bokeh.plotting import figure, curdoc
from bokeh.models.widgets import Select, Slider, RadioGroup
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application

import pypospack.visualization as visualization


'''
score += 1 for err less than mean
Histogram of scores
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

    # fix this later use length of name list
    bin_list = list(range(0, 11))
    bin_list = [v+0.5 for v in bin_list]
    hist, edges = np.histogram(total_df[name], density=True, bins=bin_list, range=(minimum, maximum))
    x = np.linspace(minimum, maximum, len(total_df[name]))
    # probability density function
    pdf = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
    pdf = pdf[np.logical_not(np.isnan(pdf))]
    gaussian_dict['hist'] = hist
    gaussian_dict['edges'] = edges
    gaussian_dict['x'] = x
    gaussian_dict['pdf'] = pdf

    return gaussian_dict


def get_score_by_param(df, name_list):

    name_avg_dict = {}
    name_score_dict = {}
    for names in name_list:
        name_avg_dict[names] = (df[names].mean())**2

    for names in name_list:
        name_score_dict[names] = 0
        for i in range(len(df)):
            if (df[names].iloc[i])**2 < name_avg_dict[names]:
                name_score_dict[names] += 1


def get_score_by_simulation(df, name_list):
    avg_val_dict = {}
    row_score_list = []
    for names in name_list:
        mean = df[names].median()
        avg_val_dict[names] = mean
    for i in range(0, len(df)):
        row = df.iloc[i]
        row_score = 0
        for names in name_list:
            value = row[names]
            if value < avg_val_dict[names]:
                row_score += 1
        row_score_list.append(row_score)

    df['score'] = row_score_list
    return row_score_list


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

        # print(len(self.param_names), len(self.err_names), len(self.sq_err_names), len(self.qoi_names))


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


    def write_text_file(self):
        scores = get_score_by_simulation(self.total_pandas_df, self.sq_err_names)
        best_simulations = []
        for i, v in enumerate(scores):
            if v == 10:
                best_simulations.append(i)
        best_sim_values = []
        for sims in best_simulations:
            pandas_row = self.total_pandas_df.iloc[sims]
            pandas_row = list(pandas_row)
            pandas_row = map(str, pandas_row)
            pandas_row = ' '.join(pandas_row)
            formatted_row = str(sims)+' '+str(pandas_row)
            best_sim_values.append(formatted_row)
        header = 'sim_id '+' '.join(list(self.total_pandas_df))
        assert len(best_sim_values) > 0
        with open('best_simulations.txt', mode='w') as f:
            f.write(header)
            f.write('\n')
            for sim_rows in best_sim_values:
                f.write(sim_rows)
                f.write('\n')


    def generate_plots(self, doc):
        self.hist_plot = {}
        get_score_by_simulation(self.total_pandas_df, self.sq_err_names)
        self.hist_plot['plotting_data'] = np.array(self.total_pandas_df['score'])
        self.hist_plot['plot_width'] = 1000
        self.hist_plot['plot_height'] = 500
        self.hist_plot['title'] = 'Simulation Score Histogram'

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
        gauss_dict = gaussian_fit(name='score', total_df=self.total_pandas_df)
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

        self.hist_plot['object_figure'].yaxis.axis_label = 'Relative Frequency'
        self.hist_plot['object_figure'].xaxis.axis_label = 'Simulation Score'

        doc.add_root(self.hist_plot['object_figure'])
        self.write_text_file()


if __name__ == "__main__":
    c_i = CreateInterface()
    c_i.start_bokeh_server()
