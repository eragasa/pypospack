# imports here

import os, time

import numpy as np
import pandas as pd

import bokeh.layouts
from bokeh.io import output_file, show, curdoc, output_notebook
from bokeh.layouts import gridplot,row,column
from bokeh.models import ColumnDataSource
from bokeh.models.callbacks import CustomJS
from bokeh.models.glyphs import Circle
from bokeh.plotting import figure, curdoc
from bokeh.client import push_session
from bokeh.models.widgets import Select,PreText
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application

class VisualizationTool(object):
    """
    Attributes:
        param_names (list of str)
        qoi_names (list of str)
        err_names (list of str)
        param_df (pandas.DataFrame)
        qoi_df (pandas.DataFrame)
        err_df (pandas.DataFrame
        total_df (pandas.DataFrame)
        source(list of dataframes)

    """
    def __init__(self):
        # initialize some attributes as None
        self.param_names = None
        self.qoi_names = None
        self.err_names = None
        self.param_df = None
        self.qoi_df = None
        self.err_df = None
        self.total_df = None

        # temp attribute for debugging
        self.pe_list = None

        # configuration variables for bokeh
        bokeh_tools = ['box_select','reset']
        self.bokeh_tools = ', '.join(bokeh_tools)

    def load_data_file(self,fname):
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

    def _start_bokeh_server(self):
        self.bokeh_app = Application(FunctionHandler(self._generate_frame))
        self.bokeh_server = Server(
            {'/':self.bokeh_app},num_procs=1
        )
        self.bokeh_server.start()

        # start io loop for bokeh_server
        self.bokeh_server.io_loop.add_callback(
            self.bokeh_server.show,'/')
        self.bokeh_server.io_loop.start()

    def _make_pandas_sources_for_bokeh(self):
        self.source_list = []   # list of pandas.DataFrame

        p_e_list = list(zip(self.param_names,self.err_names))
        for names in p_e_list:
            y_name = names[1]
            x_name = names[0]
            x_data = list(self.param_df[x_name])
            y_data = list(self.err_df[y_name])
            source_df = pd.DataFrame({x_name: x_data, y_name: y_data})
            # source_list.append(ColumnDataSource(source_df))
            self.source_list.append(source_df)

    def _make_bokeh_callbacks(self):
        self.bokeh_callbacks = {}
        for i, source in enumerate(self.source_list):
            callback_key = 'source{}'.format(i)
            self.bokeh_callbacks[callback_key] = source

    def _generate_frame(self):
        #make bokeh figure
        self.graph = {}
        self.graph['plot_width'] = 600
        self.graph['plot_height'] = 400
        self.graph['title'] = 'place holder title string'
        self.graph['data'] = self.source_list[0]

        # create bokeh.plotting.figure instance
        try:
            self.graph['figure'] = figure(
                tools = self.bokeh_tools,
                plot_width = self.graph['plot_width'],
                plot_height = self.graph['plot_height'],
                title = self.graph['title']
            )
        except:
            # debugging code
            print('Problem Setting Up Figure')
            print('plot_width:',type(self.graph['plot_width']),
                    ':',self.graph['plot_width']
            )
            print('plot_height:',type(self.graph['plot_height']),
                    ':',self.graph['plot_height']
            )
            print('title',type(self.graph['title']),
                    ':',self.graph['title'])
            raise

        # plot data
        try:
            self.graph['figure'].circle(
                x = self.graph['data'].ix[:,0],
                y = self.graph['data'].ix[:,1]
                #y = self.graph['data'].ix[:,1],
                #source = self.graph['data']
            )
        except:
            print('Problem Plotting Data')
            print('x:\n',
                    '\ttype:{}\n'.format(
                        type(self.graph['data'].ix[:,0])),
                    '\tshape:{}\n'.format(
                        self.graph['data'].ix[:,0].shape)
            )
            print('y:\n',
                    '\ttype:{}\n'.format(
                        type(self.graph['data'].ix[:,1])),
                    '\ttype:{}\n'.format(
                        self.graph['data'].ix[:,1].shape),
            )
            print('source:\n'
                    '\ttype:{}\n'.format(
                        type(self.graph['data'])),
                    '\ttype:{}\n'.format(
                        self.graph['data'].shape)
            )
            raise
 
        self.graph['figure'].show()
    def _bokeh_callback(self):
        pass

class VisualizationDemo(object):
    def __init__(self):
        bokeh_tools = ['box_select','reset']
        self.bokeh_tools = ', '.join(bokeh_tools)


    def load_data_file(self,fname):
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

    def update_data(self,param_x,param_y,err_x,err_y):
        self.total_df['param_x'] = self.total_df[param_x]
        self.total_df['param_y'] = self.total_df[param_y]
        self.total_df['err_x'] = self.total_df[err_x]
        self.total_df['err_y'] = self.total_df[err_y]

    def nix(self,val,lst):
        return [x for x in lst if x!=val]

    def setup_bokeh_frame(self):
        self.source = ColumnDataSource(
                data = dict(
                    param_x = [], 
                    param_y = [],
                    err_x = [], 
                    err_y = []
                )
        )
        self.source_static = ColumnDataSource(
                data = dict(
                    param_x = [], 
                    param_y = [],
                    err_x = [], 
                    err_y = []
                )
        )

        #self.stats = PreText(
        #        text='',
        #       width=500
        #)
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
        self.param_graph['plot_width'] = 400
        self.param_graph['plot_height'] = 400
        self.param_graph['tools'] = self.bokeh_tools
        self.param_graph['obj_figure'] = figure(
                plot_width=self.param_graph['plot_width'],
                plot_height=self.param_graph['plot_height'],
                tools=self.param_graph['tools']
        )
        self.param_graph['obj_glyph'] = Circle(
                x = 'param_x',
                y = 'param_y',
                size = 1,
        )
        self.param_graph['obj_figure'].add_glyph(
                self.source,
                self.param_graph['obj_glyph']
        )

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
        self.err_graph['plot_width'] = 400
        self.err_graph['plot_height'] = 400
        self.err_graph['tools'] = self.bokeh_tools
        self.err_graph['obj_figure'] = figure(
                plot_width=self.err_graph['plot_width'],
                plot_height=self.err_graph['plot_height'],
                tools=self.err_graph['tools']
        )
        self.err_graph['obj_glyph'] = Circle(
                x = 'err_x',
                y = 'err_y',
                size = 1,
        )
        self.err_graph['obj_figure'].add_glyph(
                self.source,
                self.err_graph['obj_glyph']
        )
        # callback functions
        def param_x_select_change(attrname,old,new):
            self.param_graph['obj_x_select'].options = \
                    self.nix(new, self.param_names)
            update()

        def param_y_select_change(attrname,old,new):
            self.param_graph['obj_y_select'].options = \
                    self.nix(new, self.param_names)
            update()

        self.param_graph['obj_x_select'].on_change(
                'value',
                param_x_select_change
        )
        self.param_graph['obj_y_select'].on_change(
                'value',
                param_y_select_change
        )

        def err_x_select_change(attrname,old,new):
            self.err_graph['obj_x_select'].options = \
                    self.nix(new, self.err_names)
            update()

        def err_y_select_change(attrname,old,new):
            self.err_graph['obj_y_select'].options = \
                    self.nix(new, self.err_names)
            update()

        self.err_graph['obj_x_select'].on_change(
                'value',
                err_x_select_change
        )
        self.err_graph['obj_y_select'].on_change(
                'value',
                err_y_select_change
        )

        def update(selected=None):
            param_name_x = self.param_graph['obj_x_select'].value
            param_name_y = self.param_graph['obj_y_select'].value
            err_name_x = self.err_graph['obj_x_select'].value
            err_name_y = self.err_graph['obj_y_select'].value

            self.update_data(
                    param_name_x, param_name_y,
                    err_name_x, err_name_y
            )
            show(layout)
        
        def selection_change(attrname,old,now):
            param_name_x = self.param_graph['obj_x_select'].value
            param_name_y = self.param_graph['obj_y_select'].value
            err_name_x = self.err_graph['obj_x_select'].value
            err_name_y = self.err_graph['obj_y_select'].value
           
        self.source.on_change('selected',selection_change)

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

        update()
        #curdoc().add_root(layout)

    def start_bokeh_server(self):
        self.bokeh_app = Application(
                FunctionHandler(self.setup_bokeh_frame))
        self.bokeh_server = Server(
            {'/':self.bokeh_app},num_procs=1
        )
        self.bokeh_server.start()

        # start io loop for bokeh_server
        self.bokeh_server.io_loop.add_callback(
            self.bokeh_server.show,'/')
        self.bokeh_server.io_loop.start()

if __name__ == "__main__":

    data_dir = 'data'
    filename = 'culled_009.out'

    #------
    vizgraph = VisualizationTool()
    vizgraph.load_data_file(
        fname=os.path.join(data_dir,filename)
    )
    
    print('checking pandas dataframes')
    for df in [self.param_df,self.qoi_df,self.err_df]:
        pass

    vizgraph._make_pandas_sources_for_bokeh()
    vizgraph._generate_frame()

    #print(vizgraph.p_e_list)
    #print(vizgraph.tools)
    #assert vizgraph.tools == 'box_select, reset'

