import os
from collections import OrderedDict

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

from pypospack.workflows import LammpsThermalExpansion

class FileSampler():

    def __init__(self,
                 configuration,
                 data,
                 structure_name,
                 structure_path,
                 workflow_type,
                 workflow_definition):
        self.initialize_configuration(configuration)
        self.initialize_data(data)
        self.structure_name = structure_name
        self.strucutre_path = structure_path
        self.workflow_type = workflow_type
        self.workflow_definition = workflow_definition
        self.potential_definition = self.configuration.potential

    def initialize_configuration(self,configuration):
        if isinstance(configuration,PyposmatConfigurationFile):
            self.configuration = configuration
        elif isinstance(configuration,str):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=configuration)
        else:
            msg = ("configuration must be a path to a configuration file or an "
                   "instance of the PyposmatConfigurationFile,")
            raise TypeError(msg)

    def initialize_data(self,data):
        if isinstance(data,PyposmatDataFile):
            self.data = data
        elif isinstance(data,str):
            self.data = PyposmatDataFile()
            self.data.read(filename=data)
        else:
            msg = ("data must be a path to a data file or an instance of "
                   "PyposmatDataFile.")
            raise TypeError(msg)

    def run(self):
        for index,row in self.data.df.iterrows():
            sim_id = row['sim_id']
            print('working on sim_id:{}'.format(sim_id))
            parameters = OrderedDict([(k,row[k]) for k in self.configuration.parameter_names])

            original_path = os.getcwd()
            os.mkdir(sim_id)
            os.chdir(sim_id)
            if workflow_type == 'lmps_thermal_expansion':
                workflow = LammpsThermalExpansion(
                        structure_name=Si_structure_definition['name'],
                        structure_path=Si_structure_definition['filename'],
                        **workflow_definition)
                workflow.create_task_configurations()
                workflow.create_tasks()
                workflow.prepare_tasks(
                        potential_definition = self.potential_definition,
                        potential_parameters = parameters)
                workflow.run()
            os.chdir(original_path)
if __name__ == "__main__":
    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()

    config_fn = os.path.join(
            pypospack_root_dir,
            'data','Si__sw__data','pareto_optimization_unconstrained',
            'pyposmat.config.in')
    data_fn = 'pyposmat.reference.in'
    
    Si_structure_definition = OrderedDict()
    Si_structure_definition['name'] = 'Si_dia_unit'
    Si_structure_definition['filename'] = os.path.join(
            os.path.join(
                pypospack.utils.get_pypospack_root_directory(),
                'data','Si__structure_db',
                'Si_dia_unit.vasp')
    )

    workflow_type = 'lmps_thermal_expansion'
    workflow_definition = OrderedDict()
    workflow_definition['temperature_low'] = 100
    workflow_definition['temperature_high'] = 1000
    workflow_definition['temperature_step'] = 100
    workflow_definition['pressure'] = 1
    workflow_definition['time_step'] = 0.001
    workflow_definition['time_total'] = 0.001 * 10000
    workflow_definition['supercell'] = [10,10,10]
    
    
    o = FileSampler(configuration=config_fn,
                    data=data_fn,
                    structure_name = Si_structure_definition['name'],
                    structure_path = Si_structure_definition['filename'],
                    workflow_type=workflow_type,
                    workflow_definition=workflow_definition)
    o.run()
