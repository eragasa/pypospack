import copy, yaml
from collections import OrderedDict
def test__read():
    _yaml_filename_in = 'pypospack.qoi.yaml'

    from pypospack.qoi import QoiDatabase
    qoidb = QoiDatabase()
    qoidb.read(filename=_yaml_filename_in)

    assert isinstance(qoidb.qois,dict)
    assert isinstance(qoidb.qois,OrderedDict)

def test__write():
    _yaml_filename_in = 'pypospack.qoi.yaml'
    _yaml_filename_out = 'pypospack.qoi.yaml.out'

    from pypospack.qoi import QoiDatabase
    qoidb = QoiDatabase()
    qoidb.read(filename=_yaml_filename_in)
    qoidb.write(filename=_yaml_filename_out)
    qoidb.read(filename=_yaml_filename_out)

if __name__ == '__main__':
    qoi_db = QoiDatabase()
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.a0', 
            qoi_type = 'a0_min_all',
            structures = ['MgO_NaCl'],
            target = 4.246)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c11', 
            qoi_type = 'c11',
            structures = ['MgO_NaCl'],
            target = 277.00)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c12', 
            qoi_type = 'c12',
            structures = ['MgO_NaCl'],
            target = 91.67)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c44', 
            qoi_type = 'c44',
            structures = ['MgO_NaCl'],
            target = 144.01)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.B', 
            qoi_type = 'bulk_modulus',
            structures = ['MgO_NaCl'],
            target = 153.45)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.G', 
            qoi_type = 'shear_modulus',
            structures = ['MgO_NaCl'],
            target = 92.66)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.fr_a', 
            qoi_type = 'defect_energy',
            structures = OrderedDict((
                ['defect','MgO_NaCl_fr_a'],
                ['reservoir','MgO_NaCl']
            )),
            target = 10.978)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.fr_c', 
            qoi_type = 'defect_energy',
            structures = OrderedDict((
                ['defect','MgO_NaCl_fr_c'],
                ['reservoir','MgO_NaCl']
            )),
            target = 8.986)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.sch', 
            qoi_type = 'defect_energy',
            structures = OrderedDict((
                ['defect','MgO_NaCl_sch'],
                ['reservoir','MgO_NaCl']
            )),
            target = 5.067)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.001s', 
            qoi_type = 'surface_energy',
            structures = ['MgO_NaCl_001_s','MgO_NaCl'],
            target = 0.05595)

    qoi_db.write(filename='pypospack.qoi.yaml')

    for k,v in qoi_db.qois.items():
        print(k,v)
