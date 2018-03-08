# imports here
'''
TODO:
scaling and box zoom

probability density plot
-histogram
1d fit
2s fit
-multivariant normal distribution

kernel density estimate
-pyflames post module

principal components analysis
clustering
manifold learning T-sne
'''
import os,copy
from pypospack.visualization import ParetoOptimizationVisualization
from pypospack.pyposmat.data import PyposmatDataFile

import numpy as np
import pandas as pd
class PypospackVisualization(ParetoOptimizationVisualization):

    def load_data_file(self, fname):
        data = PyposmatDataFile()
        data.read(filename=fname)
        self.param_names = list(data.parameter_names)
        self.qoi_names = list(data.qoi_names)
        self.err_names = list(data.error_names)
        self.names = list(
                data.parameter_names
                + data.qoi_names
                + data.error_names)
        self.df = copy.deepcopy(data.df)
        self.total_df = self.df[self.names]
        self.param_df = self.df[self.param_names]
        self.qoi_df = self.df[self.qoi_names]
        self.err_df = self.df[self.err_names]


if __name__ == "__main__":
    import argparse

    #class FullPaths(argparse.Action):
    #    """Expand user- and relative-paths"""
    #    def __call__(self, parser, namespace, values, option_string=None):
    #        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

    #def is_dir(dirname):
    #    """Checks if a path is an actual directory"""
    #    if not os.path.isdir(dirname):
    #        msg = "{0} is not a directory".format(dirname)
    #        raise argparse.ArgumentTypeError(msg)
    #    else:
    #        return dirname

    #def is_file(fn):
    #    """Checks if a path is an actual directory"""
    #    if not os.path.isfile(fn):
    #        msg = "{0} is not a directory".format(fn)
    #        raise argparse.ArgumentTypeError(msg)
    #    else:
    #        return fn

    #parser = argparse.ArgumentParser
    #parser.add_argument("--in",
    #    action=FullPaths,type=is_file,dest='filename_in',
#        help="location of the pyposmat datafile we want to analyze"
    #)
    #parser.add_argument("--out",
    #    action=FullPaths,type=is_file,dest='filename_out',
#        help="location of the subselect we want to output"
    #)
    #args = parser.parse_args()

    #fn_i = args.filename_in
    #fn_o = args.filename_out

    fn_i = "results.pareto.out"
    fn_o = "results.temp.out"
    #d = PyposmatDataFile()
    #d.read(filename=fn_i)
    #df = d.df
    #df = df.loc[
    #    (df["Ni_fcc.E_coh.err"] < 2)
    #    & (df["Ni_fcc.a0.err"] < 0.5)
    #    & (df["Ni_fcc.c11.err"] > -50)
    #    & (df["Ni_fcc.c12.err"] < 50)
    #    & (df["Ni_fcc.c12.err"] > -50)
    #    & (df["E_Ni_fcc_bcc"] > 0)
    #    & (df["E_Ni_fcc_hcp"] > 0)
    #    & (df["E_Ni_fcc_dia"] > 0)
    #    & (df["E_Ni_fcc_sc"] > 0)]
    #d.df = df
    #d.write(filename=fn_o)

    vizdemo = PypospackVisualization()
    vizdemo.load_data_file(fname=fn_i)
    vizdemo.start_bokeh_server()
    vizdemo.setup_bokeh_frame()
