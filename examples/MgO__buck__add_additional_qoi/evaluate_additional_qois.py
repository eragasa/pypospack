import os,shutil,sys,copy
import numpy as np
import pandas as pd
from mpi4py import MPI
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.engines import PyposmatIterativeSampler
from pypospack.pyposmat.engines import PyposmatEngine

class NewPyposmatFileSampler(PyposmatEngine):

    def __init__(self,
            config_fn='pyposmat.config.in',
            data_in_fn='pyposmat.results.in',
            data_out_fn='pyposmat.results.out',
            mpi_rank = None,
            mpi_size=None,
            log = None,
            base_directory = None):
        
        PyposmatEngine.__init__(self,
                filename_in=config_fn,
                filename_out=data_out_fn,
                base_directory=base_directory,
                fullauto=False)

        self.configuration_fn = config_fn
        self.datafile_in_fn = data_in_fn
        self.datafile_out_fn = data_out_fn

        self.mpi_rank = mpi_rank
        self.mpi_size = mpi_size

        self.configuration = None
        self.datafile = None
        self.subselect_df = None
        self.reference_potentials = None

        self.parameter_names = None
        self.qoi_names = None
        self.error_names = None
        self.normed_error_names = None
        
        self.qoi_validation_names = None
        self.error_validation_names = None
        self.normed_error_validation_names = None

        self.qoi_targets = None
        self.qoi_validation_target = None

        if config_fn is not None:
            self.read_configuration_file(config_fn=config_fn)

        if data_in_fn is not None:
            self.read_datafile_in(datafile_fn=data_in_fn)

        if data_out_fn is not None:
            self.datafile_out = PyposmatDataFile(data_out_fn)

    def __log(self,m):
        print(m)

    def configure_qoi_manager(self,qois=None,use_fitting_qois=True,use_testing_qois=False):

        if qois is not None:
            _qois = copy.deepcopy(qois)
        else:
            _qois = OrderedDict()

            if use_fitting_qois:
                for k,v in self.configuration.qois.items():
                    _qois[k]=v

            if use_testing_qois:
                for k,v in self.configuration.qois_validation.items():
                    _qois[k]=v
        PyposmatEngine.configure_qoi_manager(self,_qois)

    def configure_task_manager(self):
        PyposmatEngine.configure_task_manager(self)

    def read_configuration_file(self,config_fn=None):
        if config_fn is not None:
            self.configuration_fn = config_fn
        
        _config_fn = self.configuration_fn


        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=_config_fn)

        self.structure_directory = self.configuration.structures['structure_directory']
        if os.path.isdir(self.structure_directory):
            msg = "[OK] structure_directory:".format(self.structure_directory)
            self.__log(msg)
        else:
            msg = "[FAIL] structure_directory:".format(self.structure_directory)
            raise PyposmatEngineError(msg)
        
        # set parameter names
        self.parameter_names = self.configuration.parameter_names
        
        # set name arrays for testing qois
        self.qoi_names = self.configuration.qoi_names
        self.error_names = self.configuration.error_names
        self.normed_error_names = self.configuration.normed_error_names

        # set name arrays for validation qois
        self.qoi_validation_names = self.configuration.qoi_validation_names
        self.error_validation_names = self.configuration.error_validation_names
        self.normed_error_validation_names = self.configuration.normed_error_validation_names

        # set dictionaries for qoi targets
        self.qoi_targets = self.configuration.qoi_targets
        self.qoi_validation_targets = self.configuration.qoi_validation_targets
    
        # set dictionary for reference potentials
        self.reference_potentials = self.configuration.reference_potentials

    def read_datafile_in(self,datafile_fn):
        self.datafile_in = PyposmatDataFile()
        self.datafile_in.read(filename=datafile_fn)

    def configure_datafile_out(self,datafile_fn=None):
        if datafile_fn is not None:
            self.datafile_out_fn = datafile_fn

        _datafile_fn = self.datafile_out_fn

        self.datafile_out = PyposmatDataFile(_datafile_fn)

    def subselect_by_dmetric(self,nsmallest=50):

        # calculated normalized errors for qois
        for iqn,qn in enumerate(self.qoi_names):
            error_name = "{}.err".format(qn)
            normed_error_name = "{}.nerr".format(qn)

            q = self.qoi_targets[qn]
            error = self.datafile_in.df[error_name]

            self.datafile_in.df[normed_error_name] = error/q

        self.datafile_in.df['d_metric'] = np.sqrt(np.square(
            self.datafile_in.df[self.normed_error_names]).sum(axis=1))  
    
        self.subselect_df = self.datafile_in.df.nsmallest(nsmallest,'d_metric')

        return self.subselect_df
    
    def run_simulations(self):

        if self.qoip_names is not None and self.errorp_names is not None:
            self.datafile_out.write_header_section(
                filename = self.datafile_out_fn,
                parameter_names = self.parameter_names,
                qoi_names = self.qoi_names,
                error_names = self.error_names,
                qoip_names = self.qoip_names,
                errorp_names = self.errorp_names
                )
            
        if self.reference_potentials is not None:
            self._sample_from_reference_potentials(
                    reference_df = self.reference_df)

        if self.subselect_df is not None:
            self._sample_from_subselect_df(
                    subselect_df=self.subselect_df)
        else:
            self._sample_from_subselect_df(
                    subselect_df=self.datafile_in.df)

    def _sample_from_reference_potentials(self,reference_df=None):
        if reference_df is None:
            _reference_df = self.reference_df

        for idx,row in _reference_df[self.parameter_names].iterrows():
            print(idx)
            parameters = OrderedDict([(pn,row[pn]) for pn in self.parameter_names])
            results = self.evaluate_parameter_set(parameters=parameters)
            print(idx,parameters,results)
            
    def _sample_from_subselect_df(self,subselect_df=None):
        if subselect_df is None:
            _subselect_df = self.subselect_df

        for idx,row in _subselect_df[self.parameter_names].iterrows():
            parameters = OrderedDict([(pn,row[pn]) for pn in self.parameter_names])
            results = self.evaluate_parameter_set(parameters=parameters)
            #for task_name, task in self.tasks.items():
            #    print('sim_id:{}'.format(row['sim_id']))

if __name__ == "__main__":

    n_smallest = 10

    config_directory = "./data"
    config_fn = os.path.join(config_directory,'pyposmat.config.in')
    
    data_directory = "../../data/MgO_pareto_data"
    datafile_in_fn = os.path.join(data_directory,'culled_005.out')

    output_directory = "./"
    datafile_out_fn = os.path.join(output_directory,'qoiplus_005.out')

    o_config=PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_sampler = NewPyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn
            )

    print(80*'-')
    print('check reference potentials')
    print(80*'-')
    reference_potential_names = ["LC","BG1","BG2"]
    for v in reference_potential_names:
        print("{}:{}".format(
            v,
            v in o_sampler.reference_potentials
            ))
        assert v in o_sampler.reference_potentials

    o_sampler.create_base_directories()
    o_sampler.read_configuration_file()
    


    print(80*'-')
    print("[INFO] configure QoiManager")
    print(80*'-')

    def setup__configure_qoi_manager():
        config_directory = "./data"
        config_fn = os.path.join(config_directory,'pyposmat.config.in')
        
        data_directory = "../../data/MgO_pareto_data"
        datafile_in_fn = os.path.join(data_directory,'culled_005.out')

        output_directory = "./"
        datafile_out_fn = os.path.join(output_directory,'qoiplus_005.out')
        
        o_sampler = NewPyposmatFileSampler(
                config_fn=config_fn,
                data_in_fn=datafile_in_fn,
                data_out_fn=datafile_out_fn)
        o_sampler.create_base_directories()
        o_sampler.read_configuration_file()

    def test__configure_qoi_manager__from_OrderedDict():
        setup__configure_qoi_manager()

        _qois = OrderedDict()
        for k,v in o_sampler.configuration.qois_validation.items():
            _qois[k] =v

        o_sampler.configure_qoi_manager(qois=_qois)
    
        assert type(o_sampler.qoi_manager) is QoiManager
        assert type(o_sampler.qoi_manager.qoidb) is QoiDatabase
    
        for k,v in _qois.items():
            assert k in o_sampler.qoi_manager.qoi_db
   
    o_sampler.configure_qoi_manager(use_fitting_qois=False,use_testing_qois=True)
    
    #o_sampler.configure_qoi_manager(use_fitting_qois=True,use_testing_qois=False)

    from pypospack.qoi import QoiManager,QoiDatabase
    print("\to_sampler.qoi_manager:{}".format(type(o_sampler.qoi_manager))) 
    print("\to_sampler.qoi_manager.qoidb:{}".format(type(o_sampler.qoi_manager.qoidb)))
    assert type(o_sampler.qoi_manager) is QoiManager
    assert type(o_sampler.qoi_manager.qoidb) is QoiDatabase

    print()
    print('QOI DATABASE FOR THIS ITERATION')
    for qoi_name,qoi_config in o_sampler.qoi_manager.qoidb.qois.items():
        print(qoi_name)
        for k,v in qoi_config.items():
            print("\t{}:{}".format(k,v))

    print(80*'-')
    print("[INFO] configure task_manager")
    print(80*'-')

    def test__configure_task_manager(): pass

    o_sampler.configure_task_manager()
    print("\to_sampler.qoi_manager.tasks")
    for task_name,task_configuration in o_sampler.qoi_manager.tasks.items():
        print("\t\t{}".format(task_name))
        for k,v in task_configuration.items():
            print("\t\t\t{}:{}".format(k,v))

    print("\to_sampler.task_manager:{}".format(type(o_sampler.task_manager)))
    print("\to_sampler.task_manager.obj_Task:{}".format(type(o_sampler.task_manager.obj_Task)))
    for k,v in o_sampler.task_manager.obj_Task.items():
        print("\t\t{}:{}".format(k,type(v)))
    from pypospack.task import TaskManager
    assert type(o_sampler.task_manager) is TaskManager

    print(80*'-')
    print("[INFO] configure_datafile_out")
    o_sampler.configure_datafile_out()
    print(80*'-')
   
   
    print(80*'-')
    print("[INFO] subselect_by_dmetric")
    print(80*'-')
    o_sampler.subselect_by_dmetric(nsmallest=n_smallest)
    print("\to_sampler.subselect_df:{}".format(type(o_sampler.subselect_df))) 
 
    import pandas as pd
    assert type(o_sampler.subselect_df) is pd.DataFrame

    print(80*'-')
    print("[INFO] sample from subselect_df")
    print(80*'-')
    o_sampler._sample_from_subselect_df()
    
    exit()
    n_smallest = 50

    data_directory = "../../data/MgO_pareto_data"
    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    datafile_in_fn = os.path.join(data_directory,'culled_005.out')

    output_directory = "./"
    datafile_out_fn = os.path.join(output_directory,'qoiplus_005.out')

    config=PyposmatConfigurationFile()
    config.read(filename=config_fn)
    
    datafile=PyposmatDataFile()
    datafile.read(filename=datafile_in_fn)

    o_sampler = PyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn
            )
    o_sampler.create_base_directories()
    o_sampler.read_configuration_file()
    o_sampler.configure_qoi_manager()
    o_sampler.configure_task_manager()
    o_sampler.configure_pyposmat_datafile_out()

    o_sampler.subselect_by_dmetric(nsmallest=n_smallest)

    o_sampler._sample_from_subselect_df()
    for pn in ['LC','BG1','BG2']:
        o_sampler.add_reference_potential(
            potential_id=pn,
            potential_parameters=ref_parameters[pn],
            potential_qois=ref_qoi[pn]
            )

    # calculate normed errors
    qoi_names = config.qoi_names
    qoi_targets = config.qoi_targets

    for iqn,qn in enumerate(qoi_names):
        error_names = "{}.err".format(qn)
        normed_error_names = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        datafile.df[normed_error_names] = datafile.df[qn]/q-1

    normed_error_names = ['{}.nerr'.format(q) for q in qoi_names]

    # calculate distance metric for best potential selection
    datafile.df['d_metric'] = np.sqrt(np.square(datafile.df[normed_error_names]).sum(axis=1))
     
    subselect_df = datafile.df.nsmallest(n_smallest,'d_metric')
    
    # add reference data
    for pn in ['LC','BG1','BG2']:
        # add sim_id to the data row
        new_ref_data_row = OrderedDict()
        new_ref_data_row['sim_id'] = pn
        
        # add parameters to the data row
        for k,v in ref_parameters[pn].items():
            new_ref_data_row[k] = v

        for k,v in ref_qoi[pn].items():
            new_ref_data_row[k] = v
  
        # add errors to the data row
        for iqn,qn in enumerate(qoi_names):
            qhat = new_ref_data_row[qn]
            q = config.qoi_targets[qn]
            error_name = "{}.err".format(qn)
            error = qhat-q
            new_ref_data_row[error_name]=error
       
        # add the data row to the pandas dataframe
        subselect_df = subselect_df.append(
            pd.Series(new_ref_data_row),
            ignore_index=True
            )

   
