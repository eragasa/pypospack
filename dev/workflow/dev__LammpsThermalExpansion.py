import os
from collections import OrderedDict

import pypospack.utils

from lammps_thermal_expansion import LammpsThermalExpansion
# Stillinger and Weber, Phy. Rev. B. v. 31, p.5262 (1985)
Si__stillingerweber = OrderedDict()
Si__stillingerweber['potential_definition'] = OrderedDict()
Si__stillingerweber['potential_definition']['potential_type'] = 'stillingerweber'
Si__stillingerweber['potential_definition']['symbols'] = ['Si']
Si__stillingerweber['parameters'] = OrderedDict([
    ('SiSiSi_epsilon',2.1686),
    ('SiSiSi_sigma',2.0951),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',21.0),
    ('SiSiSi_gamma',1.20),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',7.049556277),
    ('SiSiSi_B',0.6022245584),
    ('SiSiSi_p',4.0),
    ('SiSiSi_q',0.0),
    ('SiSiSi_tol',0.0)
])

# the structures required for the simulation, which would come from the simulation database
Si_structure_definition = OrderedDict()
Si_structure_definition['name'] = 'Si_dia_unit'
Si_structure_definition['filename'] = os.path.join(
        os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__structure_db',
            'Si_dia_unit.vasp')
)

workflow_definition = OrderedDict()
workflow_definition['temperature_low'] = 100
workflow_definition['temperature_high'] = 1000
workflow_definition['temperature_step'] = 100
workflow_definition['pressure'] = 0
workflow_definition['time_step'] = 0.001
workflow_definition['time_total'] = 0.001 * 10000
workflow_definition['supercell'] = [10,10,10]
workflow = LammpsThermalExpansion(
        structure_name=Si_structure_definition['name'],
        workflow_path="./",
        **workflow_definition)
print(workflow.temperatures)
