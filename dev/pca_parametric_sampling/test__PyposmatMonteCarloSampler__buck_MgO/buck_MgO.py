from collections import OrderedDict
from pypospack.qoi import QoiDatabase
MgO_LewisCatlow = OrderedDict()
MgO_LewisCatlow['potential'] = OrderedDict()
MgO_LewisCatlow['potential']['potential_type'] = 'buckingham'
MgO_LewisCatlow['potential']['symbols'] = ['Mg','O']
MgO_LewisCatlow['parameters'] = OrderedDict()
MgO_LewisCatlow['parameters']['chrg_Mg'] = +2.0
MgO_LewisCatlow['parameters']['chrg_O']  = -2.0
MgO_LewisCatlow['parameters']['MgMg_A']   = 0.0 
MgO_LewisCatlow['parameters']['MgMg_rho'] = 0.5
MgO_LewisCatlow['parameters']['MgMg_C']   = 0.0
MgO_LewisCatlow['parameters']['MgO_A']    = 821.6
MgO_LewisCatlow['parameters']['MgO_rho']  = 0.3242
MgO_LewisCatlow['parameters']['MgO_C']    = 0.0
MgO_LewisCatlow['parameters']['OO_A']     = 2274.00 
MgO_LewisCatlow['parameters']['OO_rho']   = 0.1490
MgO_LewisCatlow['parameters']['OO_C']     = 27.88


MgO_qoi_db_gga = QoiDatabase()
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.246)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=277.00)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=91.67)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=144.01)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=153.45)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.G',
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=92.66)
if False:
    MgO_qoi_db_gga.add_qoi(
            qoi_name='MgO_NaCl.fr_a',
            qoi_type='point_defect',
            structures=OrderedDict([
                ('defect','MgO_NaCl_fr_a'),
                ('ideal','MgO_NaCl')]),
            target=10.978)
    MgO_qoi_db_gga.add_qoi(
            qoi_name='MgO_NaCl.fr_c',
            qoi_type='point_defect',
            structures=OrderedDict([
                ('defect','MgO_NaCl_fr_c'),
                ('ideal','MgO_NaCl')]),
            target=8.986)
    MgO_qoi_db_gga.add_qoi(
            qoi_name='MgO_NaCl.sch',
            qoi_type='point_defect',
            structures=OrderedDict([
                ('defect','MgO_NaCl_sch'),
                ('ideal','MgO_NaCl')]),
            target=5.067)
    MgO_qoi_db_gga.add_qoi(
            qoi_name='MgO_NaCl.001s',
            qoi_type='surface_energy',
            structures=OrderedDict([
                ('slab','MgO_NaCl_001s'),
                ('ideal','MgO_NaCl')]),
            target=0.05595)
