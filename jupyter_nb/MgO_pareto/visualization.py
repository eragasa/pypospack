import os
import time
import numpy as np
import pandas as pd
import sys
from bokeh.io import output_file, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
# sys.path.append(os.path.abspath(r"C:\Users\Seaton\repos\pypospack\\"))
import pypospack.potfit as potfit

def generate_frame():   # use later to optimize data frame creation
    pass

if __name__ == "__main__":
    data_dir = 'data'
    filename = 'culled_009.out'
    starttime = time.time()
    # configure the fitting engine
    #eip_config = potfit.AbstractFittingEngine(\
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
            data.append(data)

    assert len(names) == len(name_types)
    '''beginning of Seaton's work:
        most of this should be condensed to a single function
    '''
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

    # block below organizes data by type into lists
    param_data = []
    qoi_data = []
    err_data = []
    for key in param_key_index:
        param_col = data[key]
        param_data.append(param_col)
    for key in qoi_key_index:
        qoi_col = data[key]
        qoi_data.append(qoi_col)
    for key in err_key_index:
        err_col = data[key]
        err_data.append(err_col)

    print(type(param_data))
    timestart = time.time()
    param_data = np.array(param_data)
    print("time_to_cast to np.array: {}".format(\
        str(time.time()-timestart)))
    qoi_data = np.array(qoi_data)
    err_data = np.array(err_data)

    param_data = param_data[0:200,:]
    qoi_data = qoi_data[0:200,:]
    err_data = err_data[0:200,:]



    print(param_data.shape)
    print(qoi_data.shape)
    print(err_data.shape)
    exit()
    # block below converts data lists to pandas Series and declares an empty data frame for each data type
    param_sr_list = []
    qoi_sr_list = []
    err_sr_list = []



    for i, p in enumerate(param_data):
        param_sr = pd.Series(data=p[i], name=param_names[i])
        param_sr_list.append(param_sr)
    param_df = pd.DataFrame(data=None, columns=param_names)

    for i, p in enumerate(qoi_data):
        qoi_sr = pd.Series(data=p[i], name=qoi_names[i])
        qoi_sr_list.append(qoi_sr)
    qoi_df = pd.DataFrame(data=None, columns=qoi_names)

    for i, p in enumerate(err_data):
        err_sr = pd.Series(data=p[i], name=err_names[i])
        err_sr_list.append(err_sr)
    err_df = pd.DataFrame(data=None, columns=err_names)

    # block below converts series to data frames and merges them together
    df_list = []
    for sr in param_sr_list:
        df = sr.to_frame()
        df_list.append(df)
    param_df = pd.concat(df_list, axis=1)
    df_list = []
    for sr in qoi_sr_list:
        df = sr.to_frame()
        df_list.append(df)
    qoi_df = pd.concat(df_list, axis=1)
    df_list = []
    for sr in err_sr_list:
        df = sr.to_frame()
        df_list.append(df)
    err_df = pd.concat(df_list, axis=1)


    print('param df shape:', param_df.shape)
    print('qoi df shape:', qoi_df.shape)
    print('err df shape:', err_df.shape, '\n')

    output_file("visualization_test.html")
    tools = "box_select, reset"
    print(80*'-')
    start_time_ColumnDataSource = time.time()
    source = ColumnDataSource(data=[param_df, qoi_df, err_df])
    total_time_ColumnDataSource = time.time() - start_time_ColumnDataSource
    print(total_time_ColumnDataSource)


    # end of dataprocessing
    datadonetime = time.time() - starttime
    print(datadonetime)
    graph1 = figure(tools=tools, plot_width=300, plot_height=300, title=None)
    graph1.circle(param_df[param_names[0]], qoi_df[qoi_names[0]], source=source)
    graph1plottime = graph1plottime - datadonetime
    print("time to printe graph 1 is {} second".format(graph1plottime))
    graph2 = figure(tools=tools, plot_width=300, plot_height=300, title=None)
    graph2.circle(param_df[param_names[1]], err_df[err_names[0]], source=source)

    p = gridplot([[graph1, graph2]])
    show(p)

    print('names:', names)
    print('name_types:', name_types)
    print('data:\n\tnumber of entries:', len(data))

