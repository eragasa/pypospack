import copy,yaml
from collections import OrderedDict
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import QoiDatabase
from pypospack.io.filesystem import OrderedDictYAMLLoader  

import MgO

# <---------------- making a configuration file
MgO_qoi_db = QoiDatabase()
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.Ecoh',
        qoi_type='Ecoh_min_all',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.5
    )
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([
            ('ideal','MgO_NaCl')]),
        target=4.5
    )
MgO_potential = OrderedDict()
MgO_potential['potential_type'] = 'buckingham'
MgO_potential['symbols'] = ['Mg','O']
MgO_potential['cutoff_global'] = 10.0
MgO_structures = OrderedDict()
MgO_structures['structure_directory'] = 'test__PyposmatEngine'
MgO_structures['structures'] = OrderedDict()
MgO_structures['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'
MgO_configuration = PyposmatConfigurationFile()
MgO_configuration.qois = MgO_qoi_db.qois
MgO_configuration.potential = MgO_potential
MgO_configuration.structures = MgO_structures
assert isinstance(MgO_configuration.configuration,OrderedDict)
MgO_configuration.write(filename='pypospack.config.in')
MgO_configuration.read(filename='pypospack.config.in')
# <---------------- end make configuration file

engine = PyposmatEngine(
        filename_in = 'pypospack.config.in',
        filename_out = 'pypospack.config.out')
# <---------------- printout for debugging purposes
print('base_directory:{}'.format(engine.base_directory))
print('input_filename:{}'.format(engine.pyposmat_filename_in))
print('output_filename:{}'.format(engine.pyposmat_filename_out))
# <---------------- the steps of engine.configure() tested individually
#                   this is the step which configures the object from the
#                   configuration file
# engine.configure()
engine.create_base_directories()
engine.read_configuration_file()
engine.configure_qoi_manager()
engine.configure_task_manager()
# <
_parameters = MgO.MgO_LewisCatlow['parameters']
results = engine.evaluate_parameter_set(parameters=_parameters)

print(results)
