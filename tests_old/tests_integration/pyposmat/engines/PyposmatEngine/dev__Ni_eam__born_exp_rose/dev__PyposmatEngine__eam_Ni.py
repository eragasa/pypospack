import os
import shutil
import Ni__eam__born_exp_rose as configuration
from collections import OrderedDict

def cleanup_simulation_directories():
    sim_directories = [
        'Ni_fcc_vac.lmps_min_pos',
        'Ni_fcc.lmps_elastic',
        'Ni_fcc_100_unit.lmps_min_all',
        'Ni_fcc_111_s.lmps_min_pos',
        'Ni_fcc_isf.lmps_min_sf',
        'Ni_sc.lmps_min_all',
        'Ni_bcc.lmps_min_all',
        'Ni_fcc.lmps_min_all',
        'Ni_fcc_110_s.lmps_min_pos',
        'Ni_fcc_111_unit.lmps_min_all',
        'Ni_dia.lmps_min_all',
        'Ni_fcc_100_s.lmps_min_pos',
        'Ni_fcc_110_unit.lmps_min_all',
        'Ni_fcc_esf.lmps_min_sf',
        'Ni_hcp.lmps_min_all']

    for d in sim_directories:
        if os.path.isdir(d):
            shutil.rmtree(d)

#
symbols = ['Ni']

# THE LATTICE INFORMATION IS NORMALLY A QOI, BUT THE QOI'S ARE BURNED HERE
lattice_info = OrderedDict()
for s in symbols:
    lattice_info[s] = OrderedDict()

lattice_info['Ni']['lattice_type'] = 'fcc'
lattice_info['Ni']['cohesive_energy'] = -4.5
lattice_info['Ni']['bulk_modulus'] = 162      # in_GPa
lattice_info['Ni']['lattice_parameter'] = 3.52
a0 = lattice_info['Ni']['lattice_parameter']

# THIS IS COMPUTED INFORMATION AND IS ONLY TRUE FOR AN FCC LATTICE
lattice_type = lattice_info['Ni']['lattice_type']
if lattice_type == 'fcc':

    V = a0**3
    lattice_info['Ni']['equilibrium_volume_per_atom'] = V
    
    re = 1/(2**0.5)*a0
    lattice_info['Ni']['equilibrium_interatomic_distance'] = 1/(2**0.5)*a0 
potential_parameters = OrderedDict()
potential_parameters = OrderedDict()
potential_parameters['p_NiNi_phi0'] = 1.0
potential_parameters['p_NiNi_gamma'] = 2.0
potential_parameters['p_NiNi_r0'] = lattice_info['Ni']['equilibrium_interatomic_distance']
potential_parameters['d_Ni_rho0'] = 1.0
potential_parameters['d_Ni_beta'] = 4.0
potential_parameters['d_Ni_r0'] = lattice_info['Ni']['equilibrium_interatomic_distance']
potential_parameters['e_Ni_ecoh'] = lattice_info['Ni']['cohesive_energy']
potential_parameters['e_Ni_latticetype'] = lattice_info['Ni']['lattice_type']
potential_parameters['e_Ni_B'] = lattice_info['Ni']['bulk_modulus']
potential_parameters['e_Ni_a0'] = lattice_info['Ni']['lattice_parameter']

if __name__ == "__main__":
    from pypospack.pyposmat.engines import PyposmatEngine
    from pypospack.pyposmat.data import PyposmatDataFile
    from pypospack.pyposmat.data import PyposmatConfigurationFile

    filename_config = 'pypospack.config.in'
    configuration.write_configuration_file(filename=filename_config)

    engine = PyposmatEngine(
            filename_in = 'pypospack.config.in',
            filename_out = 'pypospack.config.out')
    engine.configure()

    _parameters = potential_parameters
    results = engine.evaluate_parameter_set(parameters=_parameters)
    print(results)

    print('cleaning up simulation directories...')
    cleanup_simulation_directories()

