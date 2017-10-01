import os
import pandas as pd
import sys
from bokeh.io import output_file, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import figure, curdoc
from bokeh.client import push_session
from bokeh.models.widgets import Select
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
import pypospack.potfit as potfit
from numpy import ndarray as ar
import numpy as np


def generateFrame(doc):  # use later to optimize data frame creation
    # create figure
    tools = "box_select, reset"
    source = ColumnDataSource(data=total_df)
    graph1 = figure(tools=tools, plot_width=600, plot_height=400, title='Params v. Error')
    graph1.circle(param_names[0], err_names[0], source=source)

    # create dropdown menu with callback
    options = []
    for i in range(len(param_names) - 1):
        options_str = str(param_names[i] + " & " + err_names[i])
        options.append(options_str)
    select = Select(value=options[0], title='Data Select', options=options)
    '''
        ---------------------------------------------------------------------
        ---------------------------------------------------------------------
        Very stuck on callback. Must be written in JavaScript.
        Try this for help:
        https://stackoverflow.com/questions/43231896/changing-source-of-plot-in-bokeh-with-a-callback
        ---------------------------------------------------------------------
        ---------------------------------------------------------------------
    '''
    select.callback = CustomJS(args=dict(source=source), code='''
    ''')

    # plot and show
    p = gridplot([[select], [graph1]])
    curdoc().add_root(p)
    show(p)

if __name__ == "__main__":
    data_dir = 'data'
    filename = 'culled_009.out'
    # configure the fitting engine
    # eip_config = potfit.AbstractFittingEngine(\
    #        fname_config_potential = 'pypospack.buckingham.yaml',
    #        fname_config_qoi = 'pypospack.qoi.yaml',
    #        fname_config_structures = 'pypospack.structure.yaml')

    fname = os.path.join(data_dir, filename)

    # persistent variables
    names = []
    name_types = []
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
    param_data = data[:, min(param_key_index):max(param_key_index)+1]
    qoi_data = data[:, min(qoi_key_index):max(qoi_key_index)+1]
    err_data = data[:, min(err_key_index):max(err_key_index)+1]
    print(param_data.shape)
    print(qoi_data.shape)
    print(err_data.shape)

    # generate pandas dataframes
    param_df = pd.DataFrame(data=param_data, columns=param_names)
    qoi_df = pd.DataFrame(data=qoi_data, columns=qoi_names)
    err_df = pd.DataFrame(data=err_data, columns=err_names)
    total_df = pd.concat([param_df, qoi_df, err_df], axis=1)

    # start bokeh server for callback access
    bokeh_app = Application(FunctionHandler(generateFrame))
    server = Server({'/': bokeh_app}, num_procs=1)
    server.start()
    server.io_loop.add_callback(server.show, "/")

    # print data info
    print('names:', names)
    print('name_types:', name_types)
    print('data:\n\tnumber of entries:', len(data))
    server.io_loop.start()