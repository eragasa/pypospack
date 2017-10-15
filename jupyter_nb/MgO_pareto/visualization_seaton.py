# imports here
'''
TODO:
scaling

principal components analysis 
clustering 
manifold learning T-sne
'''
import os, time

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
from bokeh.models.widgets import Select, PreText
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application



class VisualizationDemo(object):
    def __init__(self):
        bokeh_tools = ['box_select', 'reset']
        self.bokeh_tools = ', '.join(bokeh_tools)

    def load_data_file(self, fname):
        names = []
        names_types = []
        data = []

        # read in the line
        lines = None
        try:
            with open(fname, 'r') as f:
                lines = f.readlines()
        except:
            raise

        for i, line in enumerate(lines):
            line = line.strip()
            if i == 0:
                names = [v.strip() for v in line.split(',')]
            elif i == 1:
                name_types = [v.strip() for v in line.split(',')]
            else:
                data_line = [float(v.strip()) for v in line.split(',')]
                data.append(data_line)

        data = np.array(data)

        assert len(names) == len(name_types)
        # block below organizes data names by type into lists
        param_names = []
        param_key_index = []
        qoi_names = []
        qoi_key_index = []
        err_names = []
        err_key_index = []
        for i, v in enumerate(name_types):
            if v == 'param':
                param_names.append(names[i])
                param_key_index.append(i)
            elif v == 'qoi':
                qoi_names.append(names[i])
                qoi_key_index.append(i)
            elif v == 'err':
                err_names.append(names[i])
                err_key_index.append(i)

        # split array by data type
        param_data = data[:, min(param_key_index):max(param_key_index) + 1]
        qoi_data = data[:, min(qoi_key_index):max(qoi_key_index) + 1]
        err_data = data[:, min(err_key_index):max(err_key_index) + 1]

        # generate pandas dataframes
        self.param_df = pd.DataFrame(data=param_data, columns=param_names)
        self.qoi_df = pd.DataFrame(data=qoi_data, columns=qoi_names)
        self.err_df = pd.DataFrame(data=err_data, columns=err_names)
        self.total_df = pd.concat(
            [self.param_df,
             self.qoi_df,
             self.err_df], axis=1)

        # make copies to the class for persistence
        self.param_names = list(param_names)
        self.qoi_names = list(qoi_names)
        self.err_names = list(err_names)

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
        param_pane = bokeh.layouts.column(
            param_widgets,
            self.param_graph['obj_figure']
        )
        err_widgets = bokeh.layouts.row(
            self.err_graph['obj_x_select'],
            self.err_graph['obj_y_select']
        )
        err_pane = bokeh.layouts.column(
            err_widgets,
            self.err_graph['obj_figure']
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
            file_obj = open('selected_points.txt', 'w')
            selected_index_list = list(new['1d']['indices'])
            selected_rows = []
            for i in selected_index_list:
                data_row = self.total_df.iloc[i]
                selected_rows.append(data_row)
            formatted_rows = []
            for rows in selected_rows:
                param_x_row = self.param_graph['obj_x_select'].value+': '+str(rows[self.param_graph['obj_x_select'].value])
                param_y_row = self.param_graph['obj_y_select'].value+': '+str(rows[self.param_graph['obj_y_select'].value])
                err_x_row = self.err_graph['obj_x_select'].value+': '+str(rows[self.err_graph['obj_x_select'].value])
                err_y_row = self.err_graph['obj_y_select'].value+': '+str(rows[self.err_graph['obj_y_select'].value])
                formatted_rows.append(str(param_x_row)+' '+str(param_y_row)+' '+str(err_x_row)+' '+str(err_y_row))
            print(formatted_rows)

            with open('selected_points.txt', 'w') as f:
                for fr in formatted_rows:
                    f.write(fr+'\n')

        self.source.on_change('selected', source_callback)

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


if __name__ == "__main__":
    data_dir = 'data'
    filename = 'culled_009.out'

    vizdemo = VisualizationDemo()
    vizdemo.load_data_file(fname=os.path.join(data_dir, filename))
    vizdemo.start_bokeh_server()
    vizdemo.setup_bokeh_frame()
