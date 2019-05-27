import os
from collections import OrderedDict

import pypospack.utils
from pypospack.qoi import QoiDatabase

def get_testing_set_Si():
    testing_set = OrderedDict()
    # Stillinger and Weber,  Phys. Rev. B, v. 31, p. 5262, (1985)
    testing_set['parameters'] = OrderedDict([
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
    testing_set['structure_db'] = OrderedDict()
    testing_set['structure_db']['structure_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__structure_db'
    )
    testing_set['structure_db']['structures'] = OrderedDict()
    testing_set['structure_db']['structures']['Si_dia'] = 'Si_dia_unit.vasp'
    testing_set['structure_db']['structures']['Si_vac'] = 'Si_dia_333_vac.vasp'
    testing_set['qoi_db'] = QoiDatabase()
    testing_set['qoi_db'].add_qoi(
            qoi_name='Si_dia.vac',
            qoi_type='E_formation',
            structures=OrderedDict([
                ('defect','Si_vac'),
                ('ideal','Si_dia')]),
            target=3.6)
    testing_set['tasks'] = OrderedDict()
    testing_set['tasks']['Si_dia.lmps_min_all'] = OrderedDict([
            ('task_type', 'lmps_min_all'), 
            ('task_structure', 'Si_dia')
    ])
    testing_set['tasks']['Si_vac.lmps_min_pos'] = OrderedDict([
            ('task_type', 'lmps_min_pos'), 
            ('task_structure', 'Si_vac'), 
            ('bulk_structure', 'Si_dia')
    ])
    return testing_set
