import os
from collections import OrderedDict
from pypospack.pyposmat.engines import PyposmatEngine

filename_in = "pypospack.config.in"
filename_out = "pypospack.results.out"

a0=3.52
r0=1/(2**0.5)*a0
parameters = OrderedDict()
parameters['p_NiNi_phi0'] = 1.0
parameters['p_NiNi_gamma'] = 2.0
parameters['p_NiNi_r0'] = r0
parameters['d_Ni_rho0'] = 1.0
parameters['d_Ni_beta'] = 4.0
parameters['d_Ni_r0'] = r0
parameters['e_Ni_latticetype'] = 'fcc'
parameters['e_Ni_ecoh'] = -4.45
parameters['e_Ni_B']= 188.
parameters['e_Ni_a0'] = a0

if __name__ == "__main__":
    engine = PyposmatEngine(
            filename_in=filename_in,
            filename_out=filename_out)
    engine.read_configuration_file(filename=os.path.join('data','pyposmat.config.in'))
    engine.configure_qoi_manager()
    engine.configure_task_manager()
    engine.evaluate_parameter_set(parameters=parameters)

