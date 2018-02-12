# -*- coding: utf-8 -*-
"""This module provides pareto calculation functions

Eugene J. Ragasa, University of Florida, developed the original version
of this code.  Dmitriy Morozov, Lawrence Berkeley Labs, provided speed ups for the
Pareto versions of the code in Dec 2016.


"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

# <----- imports
import os,copy,sys, getopt
import numpy as np
import scipy.stats

def pareto_frontier_2d(Xs, Ys, maxX=True, maxY=True):
    '''
    Method to take two equally-sized lists and return just the elements which
    lie on the Pareto frontier, sorted into order.  Default behaviour is to 
    find the maximum for both X and Y, but the option is available to specify 
    maxX = False or maxY = False to find the minimum for either or both of the 
    parameters.
    
    original code: Jamie Bull,
    '''
    # Sort the list in either ascending or descending order of X
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    # Start the Pareto frontier with the first value in the sorted list
    p_front = [myList[0]]    
    # Loop through the sorted list
    for pair in myList[1:]:
        if maxY: 
            if pair[1] >= p_front[-1][1]: # Look for higher values of
                p_front.append(pair) # and add them to the Pareto frontier
        else:
            if pair[1] <= p_front[-1][1]: # Look for lower values of
                 p_front.append(pair) #  and add them to the Pareto frontier
    # Turn resulting pairs back into a list of Xs and Ys
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY

def dominates(p1, p2):
    """
    original code: Eugene J. Ragasa, UF
    """
    for x,y in zip(p1,p2):
        if y < x:
            return False
    return True

def pareto_bruteforce(pts, indices = None):
    """
    original code: Eugene J. Ragasa, UF
    """
    if indices is None:
        indices = list(range(len(pts)))
    result = []
    for i in indices:
        for j in indices:
            if i == j: continue
            if dominates(pts[j], pts[i]):
                break
        else:
            result.append(i)
    return result

def pareto_merge(lo, hi, i, dim):
    """
    original code: Dmitriy Morozon, LBL
    """
    if len(lo) == 0 or len(hi) == 0:
        return lo + hi

    survivors = set()
    for j in range(dim):
        if i == j: continue
        m = min(p[j] for p in lo)
        survivors.update(k for k in range(len(hi)) if hi[k][j] < m)

    return lo + [hi[k] for k in survivors]

def pareto(pts, indices = None, i = 0):
    """
    original code: Dmitriy Morozov, LBL
    """
    if indices is None:
        indices = list(range(len(pts)))
    l = len(indices)
    if l <= 1:
        return indices

    if l < 1000:
        return pareto_bruteforce(pts, indices)

    indices.sort(key = lambda x: pts[x][i])     # lazy: should use partition instead

    dim = len(pts[0])
    optimalLo = pareto(pts, indices[:l//2], (i + 1) % dim)
    optimalHi = pareto(pts, indices[l//2:], (i + 1) % dim)

    return pareto_bruteforce(pts, optimalLo + optimalHi)     # lazy: FIXME
    #return pareto_merge(optimalLo, optimalHi, i, dim)
    
class ParameterFileReader:
    """ Reads a file consisting of parameters

    This class expects a csv file the first line should have the name of all
    the parameters, separated by commas.  The second line on should have the
    value of each parameter, in the same order as the header line, also 
    separated by commas.

    Args:
        fname_param_in (str): filename of the parameter file
        
    Attributes:
        param_names     (python array) parameter names.  same index as params.
        params          (numpy.array[sim_id,param_id])
    """
    def __init__(self, fname_param_in = 'params.in'):
        self.fname_params = fname_param_in
        f = open(self.fname_params)
        self.lines = f.readlines()
        f.close()

        n_lines = len(self.lines)
        self.param_names = self.lines[0].strip().split(' ')
        self.params = []
        for idx in range(1,n_lines):
            self.params.append([float(num) for num in self.lines[idx].split()])

class SimulationResults(object):
    """this class processes the simulation results
  
    Args:
      n_simulations (int): number of simulations read from the output file
      qoi_type (str): supported qoi types are
          'abserr' - absolute error
    """

    def __init__(self):
        """default constructor"""
        self._supported_qoi_err_types = ['err','abserr', 'sqerr']

        # filenames    
        self.fname_log_file = "pyposmat.log"
        self.fname_sim_results = None
        self.fname_pareto = None
        self.fname_cull= None

        # initialize variables [ATTRIBUTES]
        self._qoi_err_type = 'abserr' # qoi error type
        self._qoi_err_req = None
        self._n_sims = None # numer of simulations
        self._names = None # names of the
        self._types = None
        self._param_names = [] # array of parameter names
        self._qoi_names = [] # array of qoi names
        self._err_names = []

        # results
        # set to None also indicates that calculations have not been done
        self._pareto_set_ids = None          # indexed with self._results
        self._results = None                 # numpy array of all simulation data
        self._pareto = None                  # numpy array of the pareto set
        self._cull = None                    # numpy array of the culled pareto set

        # we could save a lot of memory here by just storing the
        # ids of the pareto sets.
        self._pareto_id = None # todo
        self._culled_id = None # todo
        self.performance_requirements = {}

        # filename handles
        self._file_log = None
        self._log_format = None
        self._open_log_file()
 
    def __del__(self):
        self._close_log_file()

    @property
    def n_simulations(self): 
        return self._n_sims
  
    @n_simulations.setter
    def n_simulations(self, nsims):
        self._n_sims = nsims

    @property
    def qoi_ref(self):
        return self._qoi_ref

    @qoi_ref.setter
    def qoi_ref(self, dict_qoi):
        assert type(dict_qoi),dict
        self._qoi_ref = dict_qoi

    @property
    def qoi_err_req(self):
        return self._perf_req

    @property
    def names(self):
        """list of str: contains a list of all string in the column"""
        return self._names

    @property
    def types(self):
        """list of str: contains either 'param','qoi','err'"""
        return self._types

    @property
    def qoi_err_type(self):
        return self._qoi_err_type

    @qoi_err_type.setter
    def qoi_err_type(self, qe_type):
        assert type(qe_type),str

        # check to see if qoi_err_type is supported
        if qe_type not in self._supported_qoi_err_types:
            err_msg = 'unsupported qoi error type: {}'
            err_msg = err_msg.format(qe_type)
            raise ValueError(err_msg)

        self._qoi_err_type = qe_type

    @property
    def qoi_names(self):
        return self._qoi_names

    @property
    def err_names(self):
        return self._err_names

    @qoi_names.setter
    def qoi_names(self, qnames):
        self._qoi_names = qnames

    @property
    def parameter_names(self):
        return self._param_names

    @parameter_names.setter
    def parameter_names(self, pnames):
        self._param_names = pnames
    
    @property
    def results(self):
        """numpy.array: numpy array of results"""
        return self._results

    @property
    def pareto(self): 
        """numpy.array: numpy array of pareto results"""
        return self._pareto
 
    @property
    def culled(self): 
        """numpy.array: numpy array of pareto results"""
        return self._culled

    # <--- some functions here to read configuration files
    def read_configuration_files(self,
                                 pyposmat_config = 'pyposmat.config',
                                 pyposmat_potential = 'pyposmat.potential',
                                 pyposmat_qoi = 'pyposmat.qoi'):
        """" reads the pyposmat configuration files
        
        Args:
            pyposmat_config (str,optional): defaults to pyposmat.config
            pyposmat_potential (str,optional): defaults to pyposmat.potential
            pyposmat_qoi (str,optional): defaults to pyposmat.qoi
      
        """
        self._read_config_pyposmat(pyposmat_config)
        self._read_config_potential(pyposmat_potential)
        self._read_config_qoi(pyposmat_qoi)

    def _read_config_pyposmat(self,fname):
        self._config_pyposmat = pyposmat.PyPosmatConfigFile(fname,
                                                            True)

    def _read_config_potential(self,fname):
        self._config_potential = pyposmat.PotentialConfigFile(fname,
                                                              True)

    def _read_config_qoi(self,fname):
        self._config_qoi = pyposmat.QoiConfigFile(fname,
                                                  True)
        self._qoi_ref = self._config_qoi.qoi_ref_vals

    # SOME FUNCTIONS HERE TO DEAL WITH APPLICATION LOGGING.          
    def _open_log_file(self):
        if self._file_log is None:
            self._file_log = open(self.fname_log_file,'w')
  
    def _close_log_file(self):
        self._file_log.close()
  
    def _log(self,msg):
        if self._log_format is None:
            self._file_log.write(msg + "\n")
            print(msg)
        else:
            msg = self._log_format.format(msg + "\n")
            print(msg)
    
    def write_pareto_set(self,fname_out='pareto.out'):
        """Write the pareto set to file.

        This function prints the calculated pareto set to file.

        Args:
            fname_out(str) - the filename (default: pareto.out)
        """

        # create header
        str_names = ", ".join(self._names) + "\n"
        str_types = ", ".join(self._types) + "\n"

        # create body
        str_body = ""
        for sim_result in self._pareto:
            str_body += ", ".join([str(num) for num in sim_result]) + "\n"
          
        # write results
        f = open(fname_out,'w')
        f.write(str_names)
        f.write(str_types)
        f.write(str_body)
        f.close()

    def write_culled_set(self,fname_out='culled.out'):
        # create header
        str_names = ", ".join(self._names) + "\n"
        str_types = ", ".join(self._types) + "\n"

        # create body
        str_body = ""
        for sim_result in self._culled:
            str_body += ", ".join([str(num) for num in sim_result]) + "\n"
          
        # write results
        f = open(fname_out,'w')
        f.write(str_names)
        f.write(str_types)
        f.write(str_body)
        f.close()
    
    def write_analysis_files(self,
                             dir_name = None,
                             fname_pareto = 'pareto.dat',
                             fname_culled = 'culled.dat',
                             is_write_pareto = True, 
                             is_write_culled_set = True):
        """
        writes a variety of analysis files
      
        Arguments:
      
        dir_name (str) - destination directory name in which to put files
        """

        if not (dir_name == None):
            self.working_path = dir_name
        else:
            # self.working_path stays the same
            pass

        # create directory if directory does not exist
        os.makedirs(dir_name, exist_ok=True)
        msg = "working path: {}".format(self.working_path)
        self.__log(msg)

        # write results of the pareto set
        if is_write_pareto == True:
            fname = os.path.join(dir_name,fname_pareto)
            self.__log("writing pareto set to {}".format(fname))
            self.__write_pareto_set(fname)
          
        # write results of the culled pareto set
        if is_write_culled_set == True:
            fname = os.path.join(dir_name,fname_pareto)
            self.__log("writing culled pareto set to {}".format(fname))
            self.__write_culled_set()

    def __read_file(self, fname, file_type):
 
        # read file into memory
        f_in = open(fname,'r')
        lines_in = f_in.readlines()
        f_in.close()

        if file_type == 'results':
            # read header lines
            self._names = [n.strip() for n in lines_in[0].strip().split(',')]
            self._types = [t.strip() for t in lines_in[1].strip().split(',')]
        elif file_type == 'pareto':
            # check to see if pareto header line is the same as the pareto line
            if self._names == [n.strip() for n in lines_in[0].strip().split(',')]:
                pass
            else:
                if self._names is None:
                    errmsg = "The results file must be read before the pareto file"
                    raise RuntimeError(errmsg)
                else:
                    errmsg = "The pareto names header does not match results file"
                    raise RuntimeError(errmsg)

            # check to see if pareto types header line is the same as the pareto line
            if self._types == [t.strip() for t in lines_in[1].strip().split(',')]:
                pass
            else:
                if self._types is None:
                    errmsg = "The results file must be read before the pareto file"
                    raise RuntimeError(errmsg)
                else:
                    errmsg = "the pareto types header does not match results file"
                    raise RuntimeError(errmsg)

        results = []
        for i in range(2,len(lines_in)):
            result =  [v.strip() for v in lines_in[i].strip().split(',')]
            for j,t in enumerate(self._types):
                if t == "sim_id":
                    result[j] = int(float(result[j]))
                else:
                    # everything else is a float
                    result[j] = float(result[j])
            results.append(result)

        # convert into numpy file
        if file_type == 'results':
            self._param_names = [self._names[i] for i,v in enumerate(self._types) if v == 'param']
            self._qoi_names = [self._names[i] for i,v in enumerate(self._types) if v == 'qoi']
            self._err_names = [self._names[i] for i,v in enumerate(self._types) if v == 'err']
            self._results = np.array(results)
        elif file_type == 'pareto':
            self._pareto = np.array(results)
        elif file_type == 'culled':
            self._culled = np.array(results)

    def read_simulation_results(self,
                                fname_sims, 
                                fname_pareto = None, 
                                fname_cull = None):
        """
        read simulations results from a file into a memory.
      
        Args:
            fname_sims (str): the filename containing the simulation results from
                LAMMPS simulations
            fname_pareto (str): the filename containing the pareto set results
            fname_cull (str): the filename contain culled pareto set results
        """
      
        self.fname_sims = fname_sims
        self.__read_file(fname_sims, 'results')
        
        # remove rows that have NaN as result
        rows_to_remove = []
        for i in range(1,self._results.shape[0]):
            if np.isnan(self._results[i,:]).any():
                rows_to_remove.append(i)
        self._results = np.delete(self._results,rows_to_remove,axis=0)
      
        if fname_pareto is not None:
            self.fname_pareto = fname_pareto
            self.__read_file(fname_pareto, 'pareto')
          
        if fname_cull is not None:
            self.fname_cull = fname_cull
            self.__read_file(fname_cull, 'culled')
          
    def _create_dataset_for_pareto_analysis(self, err_names=None):
        """iCreates a dataset for pareto analysis

        This method creates a dataset necessary for pareto analysis

        Arguments:
        err_names (list of str): - contains the identifiers for the error
        """

        print("creating dataset for pareto analysis")
         
        if err_names is None:
            err_names = self._err_names
            
        # get indices of error names
        err_idx = [self._names.index(n) for n in err_names]

        # select the sim_id column and err_names columns
        results_err = self._results[:,[0] + err_idx]
        results_abs_err = np.abs(results_err)

        # make dataset 
        n_row, n_col = results_abs_err.shape
        self._pareto_dataset = []
        for i_row in range(n_row):
            self._pareto_dataset.append(Datapoint(i_row))
            for i_col in range(n_col):
                number = results_abs_err[i_row,i_col]
                self._pareto_dataset[i_row].addNumber(-number)
 
    def calculate_pareto_set(self):

        self._create_dataset_for_pareto_analysis(err_names=self._err_names)
        bruteforce_algo(self._pareto_dataset)

        # mark pareto set
        pareto_set = []
        pareto_set_ids = []
        for s in self._pareto_dataset:
            if s.paretoStatus == 1:
                pareto_set_ids.append(s.id)
                pareto_set.append(s.vec)
          
        #pareto_set = -np.array(pareto_set)
        self._pareto = self._results[pareto_set_ids,:]

    def calculate_parameter_estimates(self,param_list):
        params = copy.deepcopy(param_list)
        self.param_estimates = {}
        for param in params:
            self.param_estimates[param] = {}
            self.param_estimates[param]['all'] = {}
            self.param_estimates[param]['pareto'] = {}
            self.param_estimates[param]['pareto_cull'] = {}

    def calculate_qoi_estimates(self,qoi_keys):
        qois = copy.deepcopy(qoi_keys)
        set_types = ['all','pareto','pareto_cull']
        self.qoi_estimates = {}
        for qoi in qois:
            self.qoi_estimates[qoi] = {}
            for set_type in set_types:
                #TODO
                mean = 0
                std  = 1
                self.qoi_estimates[qoi][set_type] = {}
                self.qoi_estimates[qoi][set_type]['mean'] = mean
                self.qoi_estimates[qoi][set_type]['std'] = std
          
    #--------------------------------------------------------------------------
    # methods for calculating the culled pareto set
    #--------------------------------------------------------------------------
    def calculate_culled_set(self,
                             cull_type="percentile",pct=80.,
                             qoi_err_threshold = None):
        """
        Arguments:
        cull_type - supports the different culling of the pareto set by 
            different mechanisms.  The current mechanisms are 'percentile'
            and 'pct_error'
        pct - is a float variable.
        qoi_err_threshold - the error threshold acceptable.

        Notes:
            If qoi_error_threshold is set, the parameters cull_type and
            pct are ignored.

        Returns:
            Nothing

        Raises:
            RuntimeError: If any key in qoierr_threshold is not contained
                in the attribute error_names, it will check to see if
                the key value is contained in qoi_names and attempt to 
                change the key value.  If this is not successful, then
                a RuntimeError will be returned.
        """

        if qoi_err_threshold is not None:
             for k,v in qoi_err_threshold.items():
                 if k.endsin('.err'):
                     if k in self.error_names:
                         self.add_performance_constraint(k,v)
                     else:
                         raise RuntimeError('unknown performance constraint')
                 else:
                     if k in self.qoi_names:
                         new_k = "{}.err".format(k)
                         self.add_performance_constraint(new_k,v)
                     else:
                         raise RuntimeError('unknown performance constraint')
             self.apply_performance_constraints()
             return

        if cull_type == "percentile":
            self._calculate_culled_set_by_percentile(pct)
        elif cull_type == "pct_error":
            self._calculate_culled_set_by_percent_error(pct)
        else:
            raise RuntimeError("unknown cull_type")

    def _calculate_culled_set_by_percentile(self,pct_kept=80.):
        """
        Arguments:
        pct_kept (float, 10.0) - number between 1 and 100 indicating the 
        pct of simulations within the Pareto set which should be kept
                        
        Returns:
        
        a numpy array with observations indexed in rows, and parameters
        and quantities of interst indexed in columns.  The column index is the
        same as the array in "self.all_names".
        """

        # TODO:
        # A Newton-Ralphson method to get more accurate performance requirements
        #to prevent over culling of the Pareto set.
       
        if not(0 <= pct_kept <= 100.):
            errmsg = "pct_kept must be between 1 and 100, the value {} was passed."
            errmsg = errmsg.format(pct_kept)
            raise ValueError(errmsg)
        else:
            self.pct_kept = pct_kept
              
        err_keys = self._err_names
        self._perf_req = {}
        for err_key in err_keys:
            self._perf_req[err_key] = 0.
        n_sims, n_qoi = self._pareto.shape        
    
        # intialize variables
        pctl_threshold = 100        # searching for 100% within the Pareto set
                                    # to 0% in the pareto set
        is_culled = False           # intialize
        
        while not is_culled:
            rows_to_delete = []
            pctl_threshold -= 0.1
            # calculate percentile cutoffs
            for err_key in self._perf_req.keys():
                if pctl_threshold < 0:
                    errmsg = "While searching for the pctl_threshold, the \
                              percentile error dropped below zero resulting \
                              in an error."
                    raise ValueError(errmsg)
                else:
                    qoi_data = self.get_data(err_key, 'pareto','abserror')
                    cutoff = np.percentile(qoi_data,pctl_threshold)
                    self._perf_req[err_key] = cutoff

            # cull the pareto set by the performance requirements
            for err_key in self._perf_req.keys():        
                pareto = np.copy(self.pareto)
                for idx in range(n_sims):
                    ps = pareto[idx,:]

                    # determine if row needs to be deleted
                    is_delete_row = False
                    for qoi_name in self._perf_req.keys():
                        qoi_idx = self._names.index(qoi_name)
                        if ps[qoi_idx] > self._perf_req[qoi_name]:
                            is_delete_row = True 
                            
                    # add row for deletion if requirements met.
                    if is_delete_row:
                        rows_to_delete.append(idx)
            
            # check to see if the pareto set has been sufficiently culled
            n_culled = len(rows_to_delete)
            pct_culled = float(n_culled)/float(n_sims)
            if pct_kept/100. > 1 - pct_culled:
                is_culled = True

        self._culled = np.delete(self._pareto,
                                 rows_to_delete,
                                 axis=0)                
        return self._culled.copy()
            
    def _calculate_culled_set_by_percent_error(self,pct_kept=80.):
        """
        Keeps members of the pareto set, which are (1+pct) over the
        reference value
        
        Arguments:
            pct (float, 8.0) - number indicating the cutoff within which
                members of the Pareto set should be kept.
        """
        if (pct_kept) <= 0:
            errmsg = "pct_kept must be between 1 and 100, the value {} was passed."
            errmsg = errmsg.format(pct_kept)
            raise ValueError(errmsg)
        
        # calculate performance constraints.
        pct_kept = float(pct_kept)             # force casting into float
        pct_kept = pct_kept/100.
        self._perf_req = {k:pct_kept*v for k,v in self._qoi_ref.items()}
        
        self.apply_performance_constraints()

    def _calculate_culled_set_by_performance_requirements(self, perf_req=None):
        if perf_req is None:
            perf_req = self._perf_req
        else:
            assert type(perf_req), dict

        self.apply_performance_constraints()

    def add_performance_constraint(self,metric_name,metric_value):
        assert type(metric_name),str
        assert type(metric_value),float

        if self._perf_req is None:
            self._perf_req = {}
        self._perf_req[metric_name] = metric_value

    def apply_performance_constraints(self):
        # start with the full pareto set and then remove elements which do 
        # not meat the performane criteria
        n_sims, n_qoi = self._pareto.shape
        self._culled = np.copy(self._pareto)

        #determine which rows to delete
        rows_to_delete = []
        for idx in range(n_sims):
            ps = self._culled[idx,:]
            is_delete_row = False
            for qoi_name in self._perf_req.keys():
                err_name = qoi_name + ".err"
                err_idx = self._names.index(err_name)
                if np.abs(ps[err_idx]) > self._perf_req[qoi_name]:
                    is_delete_row = True 
            if is_delete_row:
                rows_to_delete.append(idx)
        
        # remove rows which do not meet performance criteria
        self._culled = np.delete(self._culled,
                                 rows_to_delete,
                                 axis=0)

    def get_data(self, name, ds_type, err_type = 'abserr'):
        """
        Arguments:
       
        name (str) - string of parameter or quantity of interest
        ds_type (str) - string of which dataset we are taking the data from.
          The ds_types which are supported are: all, pareto, pareto_culled
          
        Returns:
      
        a numpy array of the data asked for
        """
        idx  = self._names.index(name)

        # get data by dataset
        data = None # initialize
        if ds_type == 'results':
            data = self._results[:,idx]
        elif ds_type == 'pareto':
            data = self._pareto[:,idx]
        elif ds_type == 'culled':
            data = self._culled[:,idx]

        if self._types[idx] == 'err':
            # transform errors if necessary
            if err_type == 'err':
                # no transformation required
                data = self._results[:,idx]
            elif err_type == 'abserr':
                # transform for absolute errors
                data = np.abs(self._results[:,idx])
        else:
            # tranformation not necessary
            data = self._results[:,idx]

        return copy.deepcopy(data)
   
    def create_all_pareto_plots(self,qoi_list):
        for i, qoi_name_i in enumerate(qoi_list):         
            for j, qoi_name_j in enumerate(qoi_list):
                if i < j and i != j:
                    print("{} {}".format(qoi_name_i,qoi_name_j))
                    x_label = qoi_name_i
                    y_label = qoi_name_j
                    x_idx = self.all_names.index(x_label)
                    y_idx = self.all_names.index(y_label)
                    pareto_front = pareto_frontier_2d(self.np_all_sims[:,x_idx],
                                          self.np_all_sims[:,y_idx],
                                          maxX = False, maxY= False)
                    fig, ax = plt.subplots()
                    ax.scatter(self.np_all_sims[:,x_idx],
                           self.np_all_sims[:,y_idx],
                           label='dominated')
                    ax.scatter(self.pareto_set[:,self.qois.index(x_label)],
                           self.pareto_set[:,self.qois.index(y_label)],
                           label = 'pareto', color='y',)
                    ax.plot(pareto_front[0],
                            pareto_front[1],
                            color='r',
                            linewidth=2)
                    legend = ax.legend(loc="upper right")
                    plt.axis([min(pareto_front[0]),
                              max(pareto_front[0]),
                              min(pareto_front[1]),
                              max(pareto_front[1])])
                    plt.xlabel(x_label)
                    plt.ylabel(y_label)
                    plt.show()
      
class ParetoSet:
  def __init__(self):
    pass  

class Datapoint:
    """Defines a point in K-dimensional space"""
    def __init__(self,id):
        self.id = id # datapoint id (0,..N-1)
        self.vec = [] # the K-dim vector
        self.paretoStatus = -1 # -1=dont know, 1=pareto, 0=not pareto
        self.dominatedCount = 0 # number of datapoints that dominate this point
        self.dominatingSet = [] # set of vectors this one is dominating

    def addNumber(self,num):
        """Adds a number to one dimension of this datapoint"""
        self.vec.append(num)

    def addToDominatingSet(self,id2):
        """Add id of of dominating point"""
        self.dominatingSet.append(id2)

    def dominates(self,other):
        """Returns true if self[k]>=other[k] for all k and self[k]>other[k] for at least one k"""
        assert isinstance(other,Datapoint)
        gte=0 # count of self[k]>=other[k]
        gt=0 # count of self[k]>other[k]
        for k in range(len(self.vec)):
            if self.vec[k] >= other.vec[k]:
                gte+=1
                if self.vec[k] > other.vec[k]:
                    gt+=1
            
        return (gte==len(self.vec) and (gt>0))

    def __repr__(self):
        return self.vec.__repr__()+": "+str(self.paretoStatus)

def bruteforce_algo(dataset):
    num_pareto = 0
    
    # pairwise comparisons
    for n in range(len(dataset)):
        if np.mod(n,100) == 0:
          print("\t n={}".format(n))
        for m in range(len(dataset)):
            if dataset[m].dominates(dataset[n]):
                dataset[n].dominatedCount+=1
                dataset[m].addToDominatingSet(n)

    # find first pareto front
    for n in range(len(dataset)):
        if dataset[n].dominatedCount == 0:
            dataset[n].paretoStatus = 1
            num_pareto += 1
        else:
            dataset[n].paretoStatus = 0
                
def nondominated_sort(dataset):
    """Nondominated Sorting, generates ranking w/ higher number = better pareto front"""
    numPareto = 0

    # pairwise comparisons
    for n in range(len(dataset)):
        for m in range(len(dataset)):
            if dataset[m].dominates(dataset[n]):
                dataset[n].dominatedCount+=1
                dataset[m].addToDominatingSet(n)

    # find first pareto front
    front = []
    front2 = []
    tmpLevel = -10 # temporary value for Pareto level, will re-adjust later
    for n in range(len(dataset)):
        if dataset[n].dominatedCount == 0:
            dataset[n].paretoStatus = tmpLevel
            front.append(n)
            numPareto+=1

    # iteratively peel off pareto fronts
    while len(front) != 0:
        tmpLevel-=1
        for f in front:
            for s in dataset[f].dominatingSet:
                dataset[s].dominatedCount -= 1
                if dataset[s].dominatedCount == 0:
                    front2.append(s)
                    dataset[s].paretoStatus = tmpLevel
        front = front2 
        front2 = []

    # re-adjust pareto level
    for n in range(len(dataset)):
        oldLevel = dataset[n].paretoStatus
        if oldLevel != -1:
            dataset[n].paretoStatus = oldLevel-tmpLevel-1

    return numPareto


def create_dataset(raw_vectors):
    """Given a list of vectors, create list of datapoints"""
    dataset = []
    for k in range(len(raw_vectors)):
        for n,v in enumerate(raw_vectors[k]):
            if k == 0:
                dataset.append(Datapoint(n))
            dataset[n].addNumber(v)
    return dataset


def readfile(filename,multiplier=1.0):
    """Reads a vector file (objective values in one dimension)"""
    with open(filename,'r') as f:
        lines = f.readlines()
    vec = [multiplier*float(a.strip()) for a in lines]
    return vec

def make_histograms_params_combined(sim_results,params):
    #TODO: Implement this function
    raise NotImplementedError("this function has not been implemented yet")
    #print("making combined histograms for parameters")
    #n = len(params)
    #for i in range(n):
    #    param_label = params[i]
        
        
def make_histograms_qois_combined(sim_results,qois):
    print("making combined histograms")
    n = len(qois)
    for i in range(n):
        qoi_label = qois[i]
        qoi_value = sim_results.qoi_values[qoi_label]
        qoi_idx = sim_results.all_names.index(qoi_label)

        # all simulated points
        qoi_values_all = np.copy(sim_results.np_all_sims[:,qoi_idx])
        qoi_values_all = qoi_values_all[~pyflamestk.pareto.is_outlier(qoi_values_all)]
        qoi_mean_all = np.mean(qoi_values_all)
        qoi_std_all  = np.std(qoi_values_all)

        # points in the pareto set
        qoi_values_pareto = np.copy(sim_results.np_pareto_set[:,qoi_idx])
        qoi_values_pareto = qoi_values_pareto[~pyflamestk.pareto.is_outlier(qoi_values_pareto)]
        qoi_mean_pareto = np.mean(qoi_values_pareto)
        qoi_std_pareto = np.std(qoi_values_pareto)

        print("{}, mean: {}, std:{}".format(qoi_label,qoi_mean_all,qoi_std_all))
        fname = "hist_qoi_{}.jpg".format(qoi_label)
        print("filename: {}".format(fname))
        plt.figure()
        # this plot has two superimposed histograms, the first histogram
        # contains all simulated predictions, the second histogram only has
        # predictions which are part of the pareto set
        fig, ax1 = plt.subplots()
        common_params = dict(bins=20,
                             range=(min(qoi_values_all),max(qoi_values_all)),
                             alpha=1.0)
        qoi_n, qoi_bins, qoi_patches = plt.hist(qoi_values_all,    facecolor='g', **common_params)
        qoi_n, qoi_bins, qoi_patches = plt.hist(qoi_values_pareto, facecolor='y', **common_params)
        plt.xlabel(qoi_label)
        plt.axvline(x=qoi_value,color='black',linewidth=2)
        plt.ylabel('frequency')

        # add density function
        ax2 = ax1.twinx()   
        ax2.plot(qoi_bins,
                 mlab.normpdf(qoi_bins,qoi_mean_all,qoi_std_all),
                 'r--', linewidth=2)
        ax2.set_ylabel('probability density')
        plt.show()

        fname = "hist_qoi_{}.jpg".format(qoi_label)
        print("{}, mean: {}, std:{}".format(qoi_label,qoi_mean_pareto,qoi_std_pareto))
        print("filename = {}".format(fname))
        plt.figure()
        # this plot only contains one histogram contain the predictions which
        # are part of the pareto set
        fig, ax1 = plt.subplots()

        qoi_n, qoi_bins, qoi_patches = plt.hist(qoi_values_pareto, bins=20, facecolor='y', alpha=1.0)
        plt.axvline(x=qoi_value,color='black',linewidth=2)
        plt.xlabel(qoi_label)
        plt.ylabel('frequency')
        
        # add density function
        ax2 = ax1.twinx()   
        ax2.plot(qoi_bins,
                 mlab.normpdf(qoi_bins,qoi_mean_pareto,qoi_std_pareto),
                 'r--', linewidth=2)
        ax2.set_ylabel('probability density')

        plt.show()
        
def make_histograms_parameters_combined(sim_results):
    print("making parameter distribution histograms")
    
    # create dictionary item for parameter names and parameter ids
    param_names_id = {}
    for idx, param_name in enumerate(sim_results.param_names):
        if param_name != 'sim_id':
            param_names_id[param_name] = sim_results.all_names.index(param_name)
            
    # iterate over the dictionary
    for param_name, param_id in param_names_id.items():

        param_values_all = sim_results.np_all_sims[:,param_id]
        param_all_mu  = param_values_all.mean()
        param_all_std = param_values_all.std()

        param_values_pareto = sim_results.np_pareto_set[:,param_id]
        param_pareto_mu  = param_values_pareto.mean()
        param_pareto_std = param_values_pareto.std()
        
        param_values_cull = sim_results.np_pareto_set_cull[:,param_id]
        param_cull_mu  = param_values_cull.mean()
        param_cull_std = param_values_cull.std()        

        print("all:    {} {} {} {}".format(param_name, param_id, param_all_mu, param_all_std))
        print("pareto: {} {} {} {}".format(param_name, param_id, param_pareto_mu, param_pareto_std))
        print("cull:   {} {} {} {}".format(param_name, param_id, param_cull_mu, param_cull_std))

        plt.figure()
        plt.hist(param_values_all, normed=True, bins=20)
        plt.xlabel(param_name)
        plt.ylabel('probability')
        plt.title('all simulations')
        plt.show()
        
        plt.figure()
        plt.hist(param_values_pareto, normed=True, bins=20)
        plt.xlabel(param_name)
        plt.ylabel('probability')
        plt.title('pareto set')
        plt.show()

        plt.figure()
        plt.hist(param_values_cull, normed=True, bins=20)
        plt.xlabel(param_name)
        plt.ylabel('probability')
        plt.title('pareto set with performance constraints')
        plt.show()
        
def make_jpdf_plot(x_data,
                   y_data,
                   x_label,
                   y_label, 
                   axis="", 
                   title=""):
    """ Make a joint probability density plot (jpdf) plot using a kernel
    density estimate.
    
    Arguments:
    
    x_data (numpy.array)
    y_data (numpy.array)
    x_label (string)
    y_label (string)
    axis (array)
    title (string)
    """
    
    xmin = 0.
    ymax = 0.
    ymin = 0.
    ymax = 0.
    if axis == "":
        xmin = x_data.min()
        xmax = x_data.max()
        ymin = y_data.min()
        ymax = y_data.max()
        axis = [xmin,xmax,ymin,ymax]
    else:
        xmin = axis[0]
        xmax = axis[1]
        ymin = axis[2]
        ymax = axis[3]

    # prepare data for jpdf plot
    X, Y      = np.mgrid[xmin:xmax:100j,ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values    = np.vstack([x_data,y_data])
    kernel    = scipy.stats.gaussian_kde(values)
    Z         = np.reshape(kernel(positions).T, X.shape)
    
    
    plt.figure()
    plt.pcolor(X,Y,Z)
    plt.plot(x_data, y_data, 'k.', markersize=3)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.axis([xmin, xmax,ymin,ymax])
    if not title == "":
        plt.title(title)
    #plt.set_ylim([ymin, ymax])
    cb = plt.colorbar()
    cb.set_label("probability density")
    plt.show()

def make_simulations_jpdf_plot(sim_results,
                               param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                x_data  = sim_results.np_all_sims[:,x_idx]
                y_data  = sim_results.np_all_sims[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    print("Pareto plot with performance constraints")
                    print("\t{}".format(x_label,y_label))
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")    
    
def make_cull_set_jpdf_plot(sim_results,
                            param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                x_data  = sim_results.np_pareto_set_cull[:,x_idx]
                y_data  = sim_results.np_pareto_set_cull[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    print("Pareto plot with performance constraints")
                    print("\t{}".format(x_label,y_label))
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")    
                    
def make_pareto_set_jpdf_plot(sim_results,
                              param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                x_data  = sim_results.np_pareto_set[:,x_idx]
                y_data  = sim_results.np_pareto_set[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    print("Pareto plot")
                    print("\t{}".format(x_label,y_label))
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")    

def make_all_jpdf_plot(sim_results, param_list):
    for pn_i, x_label in enumerate(param_list):
        for pn_j, y_label in enumerate(param_list):
            if pn_i < pn_j and pn_i != pn_j:
                x_data  = sim_results.np_pareto_set[:,sim_results.all_names.index(x_label)]
                y_data  = sim_results.np_pareto_set[:,sim_results.all_names.index(y_label)]
                
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")
                
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                
                x_data  = sim_results.np_pareto_set_cull[:,x_idx]
                y_data  = sim_results.np_pareto_set_cull[:,y_idx]
                # make_density_plot(x_data,y_data,x_label,y_label)
                try:
                    pyflamestk.pareto.make_jpdf_plot(x_data,y_data,x_label,y_label)
                except np.linalg.linalg.LinAlgError:
                    print("Kernel Density Estimator produced LinAlgErr")

def resample_from_kernel_density(sim_results, 
                                 param_list, 
                                 param_list_not_free,
                                 fname_param = 'params.dat',
                                 n_resamples = 10000):


    # create free parameter list
    print("resampling from kernel density estimate...")
    print("\t output to file: {}".format(fname_param))
    print("\t number of resamples: {}".format(n_resamples))
    print("\t creating free parameter list...")
    free_param_list = param_list.copy()
    for param in param_list_not_free:
        free_param_list.remove(param)

    # print param list and whether or not the parameter is free        
    for param in param_list:
        if param in free_param_list:
            print("\t\t {} FREE".format(param))
        else:
            print("\t\t {} NOT FREE".format(param))

    # select only the parameters we are interested in
    param_col_index = []    
    for param_name in free_param_list:
        p_idx = sim_results.all_names.index(param_name)
        param_col_index.append(p_idx)


    #select only the parameters and do a KDE analysis
    print("\t do KDE analysis on culled parameter set...")
    cull_parameters = sim_results.np_pareto_set_cull[:,param_col_index]
    kernel = scipy.stats.gaussian_kde(cull_parameters.transpose())

    #resample from the kernel
    print("\t generation samples (may take a while)...")
    sim_results.param_sample = kernel.resample(size=n_resamples)

    #write file out
    param_list_out = []
    for idx, row in enumerate(sim_results.param_sample.T):
        row_params = []
        is_good_row = True
        for param in param_list:
            param_value = ""
            if param in free_param_list:
                param_value = row[free_param_list.index(param)]
                #print("{} in free_param_list, value = {}".format(param,param_value))
            else:
                # assignment of non-free parameters
                if param == 'chrg_O':
                    param_value = -row[free_param_list.index('chrg_Mg')]
                elif param == 'p_MgMg_a':
                    param_value = 0.
                elif param == 'p_MgMg_c':
                    param_value = 0.
                elif param == 'p_MgMg_rho':
                    param_value = 0.5
                elif param == 'p_MgO_c':
                    param_value = 0.
                else:
                    pass

            # checkparam restrictions
            if param.endswith('_rho'):
                if param_value <= 0: 
                    is_good_row = False

            row_params.append(param_value)

        if is_good_row is True:
            param_list_out.append(row_params)
    print("\t number of parameter sets: {}".format(len(param_list_out)))
    print("\t writing to file...")


    # print header
    str_out = "idx "
    for param in param_list:
      str_out += "{} ".format(param)
    
    # print params
    for idx,row in enumerate(param_list_out):
      str_out += "\n"
      str_out += "{} ".format(idx)
      for param_value in row:
          str_out += "{} ".format(param_value)
    
    
    f = open(fname_param, mode='w')
    f.write(str_out)  
    f.close()

def make_2d_pareto_plots(sim_results,
                         show_dominated=True,
                         show_pareto=True,
                         show_cull=True):
    n = len(sim_results.qoi_keys)
    for i in range(n):
        for j in range(n):
            if i < j and i != j:
                x_label = sim_results.qoi_keys[i]
                y_label = sim_results.qoi_keys[j]
                print("{} {}".format(x_label,y_label))
    
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                
                # apply performance requirements
                plt.figure()
                plt_handles = []
                if show_dominated:
                    plt_handles.append(plt.scatter(sim_results.np_all_sims[:,x_idx],
                                                   sim_results.np_all_sims[:,y_idx],
                                                   label='dominated',s=1))
                if show_pareto:
                    plt_handles.append(plt.scatter(sim_results.np_pareto_set[:,x_idx],
                                                   sim_results.np_pareto_set[:,y_idx],
                                                   color='y',
                                                   label='Pareto'))
                if show_cull:
                    plt_handles.append(plt.scatter(sim_results.np_pareto_set_cull[:,x_idx],
                                                   sim_results.np_pareto_set_cull[:,y_idx],
                                                   color='r',
                                                   label='Pareto w/ constraints'))
                                                   
                if show_pareto:
                    pareto_front = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_all_sims[:,x_idx],
                                                                        sim_results.np_all_sims[:,y_idx],
                                                                        maxX = False, maxY= False)                
                    plt.plot(pareto_front[0],
                             pareto_front[1],
                             color='y',linewidth=2)
                             
                if show_cull:
                    pareto_front_cull = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_pareto_set_cull[:,x_idx],
                                                                             sim_results.np_pareto_set_cull[:,y_idx],
                                                                             maxX = False, maxY= False)
                    plt.plot(pareto_front_cull[0],
                             pareto_front_cull[1],
                             color='r',linewidth=2)
    
                plt.axis([0,
                          np.percentile(sim_results.np_pareto_set[:,x_idx],90),
                          0,
                          np.percentile(sim_results.np_pareto_set[:,y_idx],90)])
                plt.legend(handles=plt_handles)
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                plt.show()

        
