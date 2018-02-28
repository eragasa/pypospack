from pypospack.pyposmat.data import PyposmatConfigurationFile

if __name__ == "__main__":

    import Ni__eam__morse_exp_fs_0 as config
    configuration = PyposmatConfigurationFile()
    configuration.qois = config.qoi_db.qois
    configuration.qoi_constraints = config.qoi_constraints
    configuration.structures = config.structure_db
    configuration.potential = config.potential_formalism
    configuration.sampling_type = config.sampling
    configuration.sampling_distribution = config.parameter_distribution
    configuration.sampling_constraints = config.parameter_constraints
    configuration.write(filename='pyposmat.config.in')
    configuration.read(filename='pyposmat.config.in')
