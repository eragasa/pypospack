MgO_buck_potential_definition = OrderedDict()
MgO_buck_potential_definition['potential_type'] = 'buckingham'
MgO_buck_potential_definition['symbols'] = ['Mg','O']

MgO_LC_parameters = OrderedDict()
MgO_LC_parameters['chrg_Mg'] = +2.0
MgO_LC_parameters['chrg_O']  = -2.0
MgO_LC_parameters['MgMg_A']   = 0.0 
MgO_LC_parameters['MgMg_rho'] = 0.5
MgO_LC_parameters['MgMg_C']   = 0.0
MgO_LC_parameters['MgO_A']    = 821.6
MgO_LC_parameters['MgO_rho']  = 0.3242
MgO_LC_parameters['MgO_C']    = 0.0
MgO_LC_parameters['OO_A']     = 2274.00 
MgO_LC_parameters['OO_rho']   = 0.1490
MgO_LC_parameters['OO_C']     = 27.88

MgO_structure_definition = OrderedDict()
MgO_structure_definition['name'] = 'MgO_NaCl_unit'
MgO_structure_definition['filename'] = os.path.join(
        'test_WorkflowLammpsThermalExpansion',
        'MgO_NaCl_unit.gga.relax.vasp')
MgO_structure_definition['supercell'] = [10,10,10]

qoi_configuration = OrderedDict()
qoi_configuration['qoi_name'] = 'MgO_NaCl.thermal_expansion'
qoi_configuration['temp_min'] = 100.
qoi_configuration['temp_max'] = 2000.
qoi_configuration['temp_step'] = 100.

temp = qoi_configuration

exit()

task_list = OrderedDict()

for temperature in range(
qoi_configuration['MgO_NaCl.thermal_expansion']lmps_thermal_expansion'] = 0
qoi_configuration['l
npt_task_configuration = OrderedDict()
npt_task_configuration['task'] = OrderedDict()
npt_task_configuration['task']['task_name'] = 'MgO_NaCl_unit.npt.T1000'
npt_task_configuration['task']['task_directory'] = 'MgO_NaCl_unit.npt.T1000'
npt_task_configuration['task']['dt'] = 0.001 # in picoseconds
npt_task_configuration['task']['temperature'] = 1000 # in Kelvin
npt_task_configuration[MgO_LC_configuration['task_type'] = 'lmps_npt'
MgO_LC_configuration['potential'] = MgO_buck_potential_definition
MgO_LC_configuration['parameters'] = None
MgO_LC_configuration['structure'] = MgO_structure_definition

lammps_task = LammpsNptSimulation(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename)
lammps_task.on_init(configuration)
lammps_task.on_config(configuration)
lammps_task.on_ready(configuration)

time_sleep = 0.1
print("t={}".format(temperature))
while (lammps_task.process.poll() is None):
    time.sleep(_time_sleep)
lammps_task = None

