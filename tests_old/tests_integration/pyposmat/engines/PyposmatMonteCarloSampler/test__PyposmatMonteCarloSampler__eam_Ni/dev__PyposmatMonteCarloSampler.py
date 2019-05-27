from pypospack.pyposmat.engines import PyposmatMonteCarloSampler as OldPyposmatMonteCarloSampler
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
class PyposmatMonteCarloSampler(OldPyposmatMonteCarloSampler):
    pass

def write_configuration_file(config_fn):
    """
    Args:
        config_fn(str): the name of the configuration file

    """
    from pypospack.pyposmat.data import PyposmatDataFile
    import Ni__eam__morse_exp_universal as Ni__eam
    Ni_eam_configuration = PyposmatConfigurationFile()
    Ni_eam_configuration.qois = Ni_eam.Ni_qoi_db.qois
    Ni_eam_configuration.potential = Ni_eam.Ni_eam_potential_formalism
    Ni_eam_configuration.structures = Ni_eam.Ni_structure_db
    Ni_eam_configuration.sampling_type = Ni_eam.Ni_eam_sampling
    Ni_eam_configuration.sampling_distribution =Ni_eam.Ni_eam_parameter_distribution
    Ni_eam_configuration.write(filename=config_fn)
    Ni_eam_configuration.read(filename=config_fn)
if __name__ == "__main__":

    config_fn = "pypospack.config.in"
    write_configuration_file(config_fn)




    
