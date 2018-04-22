import os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

new_data_fn = 'pyposmat.results.0.out'

# filename for lattice parameter simulations
a0_dir = os.path.join('..','MgO_condition_by_a0','data')
a0_config_fn = os.path.join(a0_dir,'pyposmat.config.in')
a0_data_fn = os.path.join(a0_dir,'pyposmat.results.0.out')

# filenamne for pressure simulations
press_dir = os.path.join('..','MgO_condition_by_pressure','data')
press_config_fn = os.path.join(press_dir,'pyposmat.config.in')
press_data_fn = os.path.join(press_dir,'pyposmat.results.0.out')


# create merged file if file does not exist
if not os.path.isfile(new_data_fn):
    print('file does not exist...')
    print('merging the files:\n\t{fn1}\n\t{fn2}'.format(
        fn1=a0_data_fn,
        fn2=press_data_fn
    ))

    # read in configuration files
    a0_config = PyposmatConfigurationFile()
    a0_config.read(a0_config_fn)
    press_config = PyposmatConfigurationFile()
    press_config.read(press_config_fn)

    # read in the datafiles
    a0_data = PyposmatDataFile()
    a0_data.read(a0_data_fn)
    press_data = PyposmatDataFile()
    press_data.read(press_data_fn)

    # doing some type checking to ensure that we read the datafile and the
    # data object *.df is actually a pandas.DataFrame
    assert type(a0_data.df) is pd.DataFrame
    assert type(press_data.df) is pd.DataFrame

    print(a0_data.df.shape)
    print(press_data.df.shape)

    assert a0_config.parameter_names == press_config.parameter_names
    parameter_names = a0_config.parameter_names
    print(a0_config.qoi_names)
    print(press_config.qoi_names)
    print(a0_config.error_names)
    print(press_config.error_names)

    new_data = PyposmatDataFile()
    new_data.write_header_section(
        parameter_names = list(parameter_names),
        qoi_names = list(a0_config.qoi_names+press_config.qoi_names),
        error_names = list(a0_config.error_names+press_config.error_names),
        filename=new_data_fn)

    for i,irow in a0_data.df[parameter_names].iterrows():
        print(i)
        for j,jrow in press_data.df[parameter_names].iterrows():


            if all([(irow[p] ==jrow[p]) for p in parameter_names]):
                print("\tmatch i({})=j({})".format(i,j))
                results = OrderedDict()

                results['parameters'] = OrderedDict([
                    (
                         p,
                         a0_data.df.iloc[i,a0_data.df.columns.get_loc(p)]
                    ) for p in parameter_names
                ])

                # add the qois to the dictionary
                results['qois'] = OrderedDict()
                for q in a0_config.qoi_names:
                    results['qois'][q] = a0_data.df.iloc[
                         i,
                         a0_data.df.columns.get_loc(q)
                    ]
                for q in press_config.qoi_names:
                    results['qois'][q] = press_data.df.iloc[
                        j,
                        press_data.df.columns.get_loc(q)
                    ]

                results['errors'] = OrderedDict()
                for e in a0_config.error_names:
                    results['errors'][e] = a0_data.df.iloc[
                         i,
                         a0_data.df.columns.get_loc(e)
                    ]
                for e in press_config.error_names:
                    results['errors'][e] = press_data.df.iloc[
                        j,
                        press_data.df.columns.get_loc(e)
                    ]

                print(results)

                # write the data results
                new_data.write_simulation_results(
                    sim_id = i,
                    results = results,
                    filename=new_data_fn
                )
                break
            else:
                pass

        #print(a0_data.df.iloc[
        #    i,
        #    [a0_data.df.columns.get_loc(c) for c in parameter_names]
        #])
        #print(press_data.df.iloc[
        #     i,
        #     [press_data.df.columns.get_loc(c) for c in parameter_names]
        #])

    #a0_col_names
    #p_col_names
else:
    print('file exists...')
