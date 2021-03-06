"""pypospack.visualization

This module contains classes for visualizatio of results done in a variety
of processed defined throughout the pypospack package library

Author: Seaton Ullberg, 2017

Version Requirements:
    bokeh > 0.12.7, to avoid tornado conflicts, https://github.com/bokeh/bokeh/issues/6152
"""

import os, time,copy

import numpy as np
import pandas as pd

import bokeh.layouts
from bokeh.io import output_file, show, curdoc, output_notebook
from bokeh.layouts import gridplot, row, column
from bokeh.models import ColumnDataSource
from bokeh.models.callbacks import CustomJS
from bokeh.models.glyphs import Circle
from bokeh.plotting import figure, curdoc
from bokeh.client import push_session
from bokeh.models.widgets import Select, PreText, TextInput
from bokeh.models import Range1d
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class PyposmatBokehVisualizer(object):
    def __init__(self):
        bokeh_tools = ['box_select', 'reset', 'box_zoom', 'pan']
        self.bokeh_tools = ', '.join(bokeh_tools)

    def read_configuration(self,filename):
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)

    def read_data(self,filename):
        self.datafile = PyposmatDataFile()
        self.datafile.read(filename=filename)

        self.parameter_names = list(self.datafile.parameter_names)
        self.qoi_names = list(self.datafile.qoi_names)
        self.error_names = list(self.datafile.error_names)

        self.param_names = list(self.datafile.parameter_names)
        self.qoi_names = list(self.datafile.qoi_names)
        self.err_names = list(self.datafile.error_names)
        print("parameter names")
        print(type(self.param_names))
        for i,v in enumerate(self.param_names):print("{:3} {:<20}".format(i,v))
        print("qoi names")
        print(type(self.qoi_names))
        for i,v in enumerate(self.qoi_names):print("{:3} {:<20}".format(i,v))
        print("error_names")
        print(type(self.err_names))
        for i,v in enumerate(self.err_names):print("{:3} {:<20}".format(i,v))
        # generate pandas dataframes
        self.param_df = copy.deepcopy(self.datafile.df[self.param_names])
        self.qoi_df = copy.deepcopy(self.datafile.df[self.qoi_names])
        self.err_df = copy.deepcopy(self.datafile.df[self.err_names])
        self.total_df = pd.concat(
            [
                self.param_df,
                self.qoi_df,
                self.err_df
            ], axis=1
        )


    def update_data(self, param_x, param_y, err_x, err_y):
        self.total_df['param_x'] = self.total_df[param_x]
        self.total_df['param_y'] = self.total_df[param_y]
        self.total_df['err_x'] = self.total_df[err_x]
        self.total_df['err_y'] = self.total_df[err_y]
        self.source.data = dict(param_x=self.total_df['param_x'],
                                param_y=self.total_df['param_y'],
                                err_x=self.total_df['err_x'],
                                err_y=self.total_df['err_y'])


    def nix(self, val, lst):
        return [x for x in lst if x != val]

    def setup_bokeh_frame(self, doc):
        self.source = ColumnDataSource(
            data=dict(
                param_x=[],
                param_y=[],
                err_x=[],
                err_y=[]
            )
        )
        self.source_static = ColumnDataSource(
            data=dict(
                param_x=[],
                param_y=[],
                err_x=[],
                err_y=[]
            )
        )
        '''
        ---------------------------------------------------------------
        Define Param Graph
        ---------------------------------------------------------------
        '''
        self.param_graph = {}
        self.param_graph['obj_x_select'] = Select(
            value=self.param_names[0],
            options=self.nix(
                self.param_names[1],
                self.param_names
            )
        )

        self.param_graph['obj_y_select'] = Select(
            value=self.param_names[1],
            options=self.nix(
                self.param_names[0],
                self.param_names
            )
        )

        self.param_graph['x_min_entry'] = TextInput(placeholder='Min X Value',
                                                    value='')
        self.param_graph['x_max_entry'] = TextInput(placeholder='Max X Value',
                                                    value='')
        self.param_graph['y_min_entry'] = TextInput(placeholder='Min Y Value',
                                                    value='')
        self.param_graph['y_max_entry'] = TextInput(placeholder='Max Y Value',
                                                    value='')

        self.param_graph['plot_width'] = 610
        self.param_graph['plot_height'] = 400
        self.param_graph['tools'] = self.bokeh_tools
        self.param_graph['obj_figure'] = figure(
            plot_width=self.param_graph['plot_width'],
            plot_height=self.param_graph['plot_height'],
            tools=self.param_graph['tools'],
            title=self.param_graph['obj_x_select'].value+' vs. '+self.param_graph['obj_y_select'].value
        )
        self.param_graph['obj_figure'].xaxis.axis_label = self.param_graph['obj_x_select'].value
        self.param_graph['obj_figure'].yaxis.axis_label = self.param_graph['obj_y_select'].value
        self.param_graph['obj_glyph'] = Circle(
            x='param_x',
            y='param_y',
            size=1,
            fill_color='#5F77D5',
            line_color='#5F77D5'
        )
        self.param_graph['obj_figure'].add_glyph(
            self.source,
            self.param_graph['obj_glyph']
        )
        '''
        ---------------------------------------------------------------
        Define Err Graph
        ---------------------------------------------------------------
        '''
        self.err_graph = {}
        self.err_graph['obj_x_select'] = Select(
            value=self.err_names[0],
            options=self.nix(
                self.err_names[1],
                self.err_names
            )
        )
        self.err_graph['obj_y_select'] = Select(
            value=self.err_names[1],
            options=self.nix(
                self.err_names[0],
                self.err_names
            )
        )

        self.err_graph['x_min_entry'] = TextInput(placeholder='Min X Value',
                                                  value='')
        self.err_graph['x_max_entry'] = TextInput(placeholder='Max X Value',
                                                  value='')
        self.err_graph['y_min_entry'] = TextInput(placeholder='Min Y Value',
                                                  value='')
        self.err_graph['y_max_entry'] = TextInput(placeholder='Max Y Value',
                                                  value='')

        self.err_graph['plot_width'] = 610
        self.err_graph['plot_height'] = 400
        self.err_graph['tools'] = self.bokeh_tools
        self.err_graph['obj_figure'] = figure(
            plot_width=self.err_graph['plot_width'],
            plot_height=self.err_graph['plot_height'],
            tools=self.err_graph['tools'],
            title=self.err_graph['obj_x_select'].value + ' vs. ' + self.err_graph['obj_y_select'].value
        )
        self.err_graph['obj_figure'].xaxis.axis_label = self.err_graph['obj_x_select'].value
        self.err_graph['obj_figure'].yaxis.axis_label = self.err_graph['obj_y_select'].value
        self.err_graph['obj_glyph'] = Circle(
            x='err_x',
            y='err_y',
            size=1,
            fill_color='#5F77D5',
            line_color='#5F77D5'
        )
        self.err_graph['obj_figure'].add_glyph(
            self.source,
            self.err_graph['obj_glyph']
        )

        def update():
            param_name_x = self.param_graph['obj_x_select'].value
            param_name_y = self.param_graph['obj_y_select'].value
            err_name_x = self.err_graph['obj_x_select'].value
            err_name_y = self.err_graph['obj_y_select'].value

            self.update_data(
                param_name_x, param_name_y,
                err_name_x, err_name_y
            )

        param_widgets = bokeh.layouts.row(
            self.param_graph['obj_x_select'],
            self.param_graph['obj_y_select']
        )
        param_x_entry = bokeh.layouts.row(
            self.param_graph['x_min_entry'],
            self.param_graph['x_max_entry']
        )
        param_y_entry = bokeh.layouts.row(
            self.param_graph['y_min_entry'],
            self.param_graph['y_max_entry']
        )
        param_pane = bokeh.layouts.column(
            param_widgets,
            self.param_graph['obj_figure'],
            param_x_entry,
            param_y_entry
        )
        err_widgets = bokeh.layouts.row(
            self.err_graph['obj_x_select'],
            self.err_graph['obj_y_select']
        )
        err_x_entry = bokeh.layouts.row(
            self.err_graph['x_min_entry'],
            self.err_graph['x_max_entry']
        )
        err_y_entry = bokeh.layouts.row(
            self.err_graph['y_min_entry'],
            self.err_graph['y_max_entry']
        )
        err_pane = bokeh.layouts.column(
            err_widgets,
            self.err_graph['obj_figure'],
            err_x_entry,
            err_y_entry
        )
        layout = bokeh.layouts.row(
            param_pane,
            err_pane
        )
        doc.add_root(layout)
        update()

        # callback functions
        def param_x_select_change(attrname, old, new):
            self.source.data['param_x'] = self.total_df[new]
            self.param_graph['obj_figure'].title.text = new+' vs. '+self.param_graph['obj_y_select'].value
            self.param_graph['obj_figure'].xaxis.axis_label = new

        def param_y_select_change(attrname, old, new):
            self.source.data['param_y'] = self.total_df[new]
            self.param_graph['obj_figure'].title.text = self.param_graph['obj_x_select'].value+' vs. '+new
            self.param_graph['obj_figure'].yaxis.axis_label = new

        self.param_graph['obj_x_select'].on_change('value', param_x_select_change)
        self.param_graph['obj_y_select'].on_change('value', param_y_select_change)

        def err_x_select_change(attrname, old, new):
            self.source.data['err_x'] = self.total_df[new]
            self.err_graph['obj_figure'].title.text = new + ' vs. ' + self.err_graph['obj_y_select'].value
            self.err_graph['obj_figure'].xaxis.axis_label = new

        def err_y_select_change(attrname, old, new):
            self.source.data['err_y'] = self.total_df[new]
            self.err_graph['obj_figure'].title.text = self.err_graph['obj_x_select'].value+' vs. '+new
            self.err_graph['obj_figure'].yaxis.axis_label = new


        self.err_graph['obj_x_select'].on_change('value', err_x_select_change)
        self.err_graph['obj_y_select'].on_change('value', err_y_select_change)

        def source_callback(attrname, old, new):
            selected_index_list = list(new['1d']['indices'])
            selected_rows = []
            for i in selected_index_list:
                data_row = self.total_df.ix[i]
                selected_rows.append(data_row)
            formatted_rows = []
            for rows in selected_rows:
                rows = rows[:-4]    # remove the 4 copied columns used in source callback
                formatted_rows.append(list(rows.get_values()))
            '''
            for rows in selected_rows:
                param_x_row = self.param_graph['obj_x_select'].value+': '+str(rows[self.param_graph['obj_x_select'].value])
                param_y_row = self.param_graph['obj_y_select'].value+': '+str(rows[self.param_graph['obj_y_select'].value])
                err_x_row = self.err_graph['obj_x_select'].value+': '+str(rows[self.err_graph['obj_x_select'].value])
                err_y_row = self.err_graph['obj_y_select'].value+': '+str(rows[self.err_graph['obj_y_select'].value])
                formatted_rows.append(str(param_x_row)+' '+str(param_y_row)+' '+str(err_x_row)+' '+str(err_y_row))
            '''
            with open('selected_points.txt', 'w') as f:
                f.write(' '.join(self.param_names)+' '+' '.join(self.err_names)+'\n')
                for fr in formatted_rows:
                    # apparently python cannot write a list to a file so the extra formatting is necessary
                    fr = str(fr)
                    fr.replace('[', '')
                    fr.replace(']', '')
                    f.write(fr+'\n')
        self.source.on_change('selected', source_callback)

        def param_x_min_callback(attrname, old, new):
            self.param_graph['obj_figure'].x_range.start = float(new)

        def param_x_max_callback(attrname, old, new):
            self.param_graph['obj_figure'].x_range.end = float(new)

        def param_y_min_callback(attrname, old, new):
            self.param_graph['obj_figure'].y_range.start = float(new)

        def param_y_max_callback(attrname, old, new):
            self.param_graph['obj_figure'].y_range.end = float(new)

        self.param_graph['x_min_entry'].on_change('value', param_x_min_callback)
        self.param_graph['x_max_entry'].on_change('value', param_x_max_callback)
        self.param_graph['y_min_entry'].on_change('value', param_y_min_callback)
        self.param_graph['y_max_entry'].on_change('value', param_y_max_callback)

        def err_x_min_callback(attrname, old, new):
            self.err_graph['obj_figure'].x_range.start = float(new)

        def err_x_max_callback(attrname, old, new):
            self.err_graph['obj_figure'].x_range.end = float(new)

        def err_y_min_callback(attrname, old, new):
            self.err_graph['obj_figure'].y_range.start = float(new)

        def err_y_max_callback(attrname, old, new):
            self.err_graph['obj_figure'].y_range.end = float(new)

        self.err_graph['x_min_entry'].on_change('value', err_x_min_callback)
        self.err_graph['x_max_entry'].on_change('value', err_x_max_callback)
        self.err_graph['y_min_entry'].on_change('value', err_y_min_callback)
        self.err_graph['y_max_entry'].on_change('value', err_y_max_callback)

    def start_bokeh_server(self):
        self.bokeh_app = Application(
            FunctionHandler(self.setup_bokeh_frame))
        self.bokeh_server = Server(
            {'/': self.bokeh_app}, num_procs=1
        )
        self.bokeh_server.start()

        # start io loop for bokeh_server
        self.bokeh_server.io_loop.add_callback(
            self.bokeh_server.show, '/'
        )
        self.bokeh_server.io_loop.start()

def start_bokeh_visualization(config_fn,data_fn):
    o = PyposmatBokehVisualizer()
    o.read_configuration(filename=config_fn)
    o.read_data(filename=data_fn)
    o.start_bokeh_server()
    o.setup_bokeh_frame()

