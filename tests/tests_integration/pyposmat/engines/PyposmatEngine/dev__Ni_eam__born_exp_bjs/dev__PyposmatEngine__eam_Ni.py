import os
import shutil
import Ni__eam__born_exp_bjs as configuration

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

    # _parameters = potential_parameters
    # results = engine.evaluate_parameter_set(parameters=_parameters)
    # print(results)

    print('cleaning up simulation directories...')
    cleanup_simulation_directories()

