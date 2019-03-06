import os, sys, shutil, copy, time
import pyflamestk.pyposmat as pyposmat
import pyflamestk.pareto as pareto
import pyflamestk.paretopost as paretopost
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.tools.plotting import scatter_matrix
import scipy.stats

class ParetoFittingPostProcessor(object):
    """
    Attributes:
        n_iterations (int)
        free_param_names (:obj:`list` of :obj:`str`)
        qoi_ref_values (:obj:`list` of :obj:`float`)
        data_dir (str)
        simulation_results (:obj:`list` of :obj:`pyflamestk.pareto.SimulationResults`)
        err_idx (:obj:`list` of :obj:`int`)
        qoi_idx (:obj:`list` of :obj:`int`)
        param_idx (:obj:`list` of :obj:`int`)
        free_param_idx (:obj:`list` of :obj:`int)

    """
    filename_results_format = "results_{:03d}.out"
    filename_pareto_format = "pareto_{:03d}.out"
    filename_culled_format = "culled_{:03d}.out"

    def __init__(self,n_iterations,data_dir):
        self.n_iterations = n_iterations
        self.free_param_names = None
        self.qoi_ref_values = None
        self.data_dir = data_dir
        self.simulation_results = []

    def load_simulation_results(self):
        self.sim_results = []
        for i in range(self.n_iterations):
            # filenames
            fname_results_in = os.path.join(\
                    self.data_dir,
                    self.filename_results_format.format(i))
            fname_pareto_in = os.path.join(\
                    self.data_dir,
                    self.filename_pareto_format.format(i))
            fname_culled_in = os.path.join(\
                    self.data_dir,
                    self.filename_culled_format.format(i))
            # load results
            self.sim_results.append(pareto.SimulationResults())
            self.sim_results[i].read_simulation_results(fname_results_in,
                                                   fname_pareto_in,
                                                   fname_culled_in)
        self.err_idx = [i for i,v in enumerate(self.sim_results[0].types) if v == 'err']
        self.qoi_idx = [i for i,v in enumerate(self.sim_results[0].types) if v == 'qoi']
        self.param_idx = [i for i,v, in enumerate(self.sim_results[0].types) if v == 'param']
        self.free_param_idx = [self.sim_results[0].names.index(n) for n in self.free_param_names]

        self.results = [[] for i in range(n_iterations)]
        self.pareto = [[] for i in range(n_iterations)]
        self.culled = [[] for i in range(n_iterations)]
        for n in range(n_iterations):
            self.results[n] = {}
            self.pareto[n] = {}
            self.culled[n] = {}

            self.results[n]['param'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].results[:,i] for i in self.free_param_idx})
            self.results[n]['qoi'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].results[:,i]\
                    for i in self.qoi_idx})
            self.results[n]['err'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].results[:,i]\
                    for i in self.err_idx})
            self.results[n]['abs_err'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:np.abs(self.sim_results[n].results[:,i])\
                    for i in self.err_idx})
            self.pareto[n]['param'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].pareto[:,i]\
                    for i in self.free_param_idx})
            self.pareto[n]['qoi'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].pareto[:,i]\
                    for i in self.qoi_idx})
            self.pareto[n]['err'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].pareto[:,i]\
                    for i in self.err_idx})
            self.pareto[n]['abs_err'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:np.abs(self.sim_results[n].pareto[:,i])\
                    for i in self.err_idx})
            self.culled[n]['param'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].culled[:,i]\
                    for i in self.free_param_idx})
            self.culled[n]['qoi'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].culled[:,i]\
                    for i in self.qoi_idx})
            self.culled[n]['err'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:self.sim_results[n].culled[:,i]\
                    for i in self.err_idx})
            self.culled[n]['abs_err'] = pd.DataFrame(\
                    {self.sim_results[n].names[i]:np.abs(self.sim_results[n].culled[:,i])\
                    for i in self.err_idx})

    def make_kde_1d_plots(self,i_iteration,subset_type,data_type):
        data = getattr(self,subset_type)[i_iteration][data_type]

