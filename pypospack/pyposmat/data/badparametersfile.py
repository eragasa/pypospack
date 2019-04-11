import os,copy
from collections import OrderedDict
import pandas as pd
import numpy as np

from pypospack.exceptions import BaseException
from pypospack.pyposmat.data import PyposmatConfigurationFile

class PyposmatBadParametersFile(object):

    def __init__(self,
                 filename = 'pyposmat.badparameters.out',
                 config_fn=None,
                 o_config=None,):
        assert filename is None or isinstance(filename,str)
        assert config_fn is None or isinstance(config_fn,str)
        assert o_config is None or isinstance(o_config,PyposmatConfigurationFile)

        self.filename = filename
        self._parameter_names = None 
        self.w_cluster_id = False

        self.df = None
        self.parameter_df = None

        self.initialize_configuration(config_fn=config_fn,o_config=o_config)

    def initialize_configuration(self,config_fn,o_config):
        if isinstance(config_fn,str) and o_config is None:
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=config_fn)
        elif isinstance(o_config,PyposmatConfigurationFile) and config_fn is None:
            self.configuration = o_config
        elif config_fn is None and o_config is None:
            self.configuration = None
        elif isinstance(o_config,PyposmatConfigurationFile) and isinstance(config_fn,str):
            m = (
                    "Cannot configure the PyposmatDataAnalyzer with both options "
                    "o_config and config_fn both being specified.  Choose only one."
            )
            raise TypeError(m)
        else:
            m = (
                    "wrong types:\n"
                    "\to_config:{}\n"
                    "\tconfig_fn:{}\n"
            )
            m = m.format(type(o_config),type(config_fn))
            raise TypeError(m)

    @property
    def parameter_names(self):

        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.parameter_names
        else:
            return None

    @property
    def names(self):

        names = ['sim_id']
        if self.df is not None:
            if 'cluster_id' in self.df.columns:
                names += ['cluster_id']
        names += self.parameter_names
        names += ['reason']

        return names

    @property
    def types(self):
        types = ['sim_id']
        if self.df is not None:
            if 'cluster_id' in self.df.columns:
                names += ['cluster_id']
        types += len(self.parameter_names)*['param']
        types += ['reason']
        
        return types

    @property
    def n_samples(self):
        (n_rows,n_cols) = self.df.shape
        return n_rows

    def get_header_string(self):

        header_line_1 = ['sim_id']\
                + self.parameter_names\
                + ['reason']
        header_line_2 = ['sim_id']\
                + len(self.parameter_names)*['param']\
                + ['reason']

        str_header_section  = "{}\n".format(",".join(header_line_1))
        str_header_section += "{}\n".format(",".join(header_line_2))

        return str_header_section

    def write_header_section(self,filename=None):

        assert isinstance(self.parameter_names,list)
        assert filename is None or isinstance(filename,str)

        if filename is not None:
            self.filename=filename

        header_str = self.get_header_string()

        with open(self.filename,'w') as f:
            f.write(header_str)

    def write_simulation_exception(self,sim_id, exception):
        is_debug = False
        assert isinstance(exception,BaseException)
        assert isinstance(self.parameter_names,list)

        if not isinstance(sim_id,str):
            sim_id=str(sim_id)
        s_reason = exception.explain()
        
        s = ",".join(
                [sim_id] \
                + [str(exception.kwargs['parameters'][k]) for k in self.parameter_names] \
                + [s_reason]
        ) + "\n"

        if is_debug:
            print(s)

        with open(self.filename,'a') as f:
            f.write(s)

    def read(self,filename=None):
        if filename is not None:
            self.filename = filename
        try:
            with open(self.filename,'r') as f:
                lines = f.readlines()
        except FileNotFoundError as e:
            print("cannot find file: {}".format(self.filename))
            print("current_working_dir: {}".format(os.getcwd()))
            raise

        self._names = [s.strip() for s in lines[0].strip().split(',')]
        self._types = [s.strip() for s in lines[1].strip().split(',')]

        table = []
        for line in lines[2:]:
            tokens = line.strip().split(',')
            values = []
            for i,v in enumerate(self._names):
                try:
                    values.append(float(tokens[i]))
                except ValueError as e:
                    if v.endswith('latticetype'):
                        values.append(tokens[i])
                    elif v.endswith('sim_id'):
                        values.append(tokens[i])
                    elif v.endswith('reason'):
                        values.append(tokens[i])
                    else:
                        print(v)
                        print(tokens[i])
                        raise
            table.append(values)

        self.df = pd.DataFrame(data=table, columns=self._names)

    def write(self,filename):
        fn = filename

        s =  [",".join(self.names)]
        s += [",".join(self.types)]
        s += [",".join([str(v) for v in k]) for k in self.df[self.names].values.tolist()]

        with open (fn, 'w') as f:
            f.write("\n".join(s))

    def write_subselect(self,filename=None):
        _filename = None
        if filename is None:
            _filename = "subselect.{}.out".format(
                    ".".join(self.score_names))
        elif type(filename) is str:
            _filename = filename
        else:
            err_msg = "the filename argument for this method must be a string"
            raise ValueError(err_msg)

        if type(self.sub_df) is None:
            err_msg="no subselection has been done on the data"
            raise ValueError(err_msg)
        if not isinstance(self.sub_df,pd.DataFrame) :
            print(type(self.sub_df))
            err_msg = "the sub_df attribute must be a pandas.DataFrame"
            raise ValueError(err_msg)

        # build the string
        str_out  = ','.join([n for n in self.names]) + "\n"
        str_out += ','.join([t for t in self.types]) + "\n"
        for row in self.sub_df.iterrows():
            _row = [a for i,a in enumerate(row[1])] #unpack tuple
            try:
                _row[0] = int(_row[0]) # row[0] is the sim_id
            except ValueError as e:
                raise
            #pass
            str_out += ','.join([str(s) for s in _row]) + "\n"

        with open(_filename,'w') as f:
            f.write(str_out)

        return _filename
