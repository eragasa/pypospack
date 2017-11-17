import os, copy
from collections import OrderedDict
import numpy as np
import pandas as pd
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatFileSampler
data_directory = os.path.join(
        '../../../test_data/test_PyposmatData',
        'data',
        'output')
pypospack_data_filename = os.path.join(
        data_directory,
        'culled_009.out')
names = [
        'chrg_Mg','chrg_O',
        'MgMg_A','MgMg_rho','MgMg_C',
        'MgO_A','MgO_rho','MgO_C',
        'OO_A','OO_rho','OO_C',
        'MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
        'MgO_NaCl.B','MgO_NaCl.G',
        'MgO_NaCl.fr_a','MgO_NaCl.fr_c','MgO_NaCl.sch',
        'MgO_NaCl.001s',
        'MgO_NaCl.a0.err','MgO_NaCl.c11.err','MgO_NaCl.c12.err',
        'MgO_NaCl.c44.err','MgO_NaCl.B.err','MgO_NaCl.G.err',
        'MgO_NaCl.fr_a.err','MgO_NaCl.fr_c.err','MgO_NaCl.sch.err',
        'MgO_NaCl.001s.err']
latex_labels = {
        'MgO_NaCl.a0':r'$a_0$',
        'MgO_NaCl.c11':r'$c_{11}$',
        'MgO_NaCl.c12':r'$c_{12}$',
        'MgO_NaCl.c44':r'$c_{44}$',
        'MgO_NaCl.B':r'$B$',
        'MgO_NaCl.G':r'$G$',
        'MgO_NaCl.fr_a':r'$E_{fr,a}$',
        'MgO_NaCl.fr_c':r'$E_{fr,c}$',
        'MgO_NaCl.sch':r'$E_{sch}$',
        'MgO_NaCl.001s':r'$\gamma_{001}$'}
parameter_names = ['chrg_Mg','chrg_O',
        'MgMg_A','MgMg_rho','MgMg_C',
        'MgO_A','MgO_rho','MgO_C',
        'OO_A','OO_rho','OO_C']
qoi_names = [
        'MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
        'MgO_NaCl.B','MgO_NaCl.G', 'MgO_NaCl.fr_a','MgO_NaCl.fr_c',
        'MgO_NaCl.sch','MgO_NaCl.001s']
error_names = ['{}.err'.format(s) for s in qoi_names]
qoi_reference_dft = {
        'MgO_NaCl.a0': 4.246,
        'MgO_NaCl.c11': 277.00031,
        'MgO_NaCl.c12': 91.67016,
        'MgO_NaCl.c44': 144.00722,
        'MgO_NaCl.B': 153.4468767,
        'MgO_NaCl.G': 92.665075,
        'MgO_NaCl.fr_a': 10.9781666,
        'MgO_NaCl.fr_c': 8.98642095,
        'MgO_NaCl.sch':5.067179685,
        'MgO_NaCl.001s': 0.055950069}
reference_LC = [2.0,-2.0,0.0,0.5,0.0,821.6,0.3242,0.0,22764.0,0.149,27.88,4.21078276561128,307.5718095128,171.135602331774,168.168424521017,216.61433805878266,68.21810359051298,9.68024989411606,9.810715180656189,5.797051474639375,0.06783817649966246,-0.035217234388720264,30.571809512799973,79.46560233177401,24.158424521016997,63.164338058782675,-24.441896409487015,-1.2977501058839405,0.8247151806561881,0.7300514746393745,0.011888176499662464]
reference_BG1 = [2.0,-2.0,0.0,0.5,0.0,1279.69,0.29969,0.0,9547.96,0.21916,32.0,4.20923604431415,383.274119165401,169.434215310753,179.601185701851,240.71418326230233,106.91995192732399,12.419259511088967,11.869114175328832,7.198887069605007,0.08070791160146304,-0.036763955685850114,106.27411916540098,77.76421531075299,35.591185701851,87.26418326230234,14.259951927323996,1.4412595110889672,2.8831141753288314,2.131887069605007,0.02475791160146304]
reference_BG2 = [1.7,-1.7,0.0,0.5,0.0,929.69,0.29909,0.0,4870,0.2679,77.0,4.222448,301.315822490901,150.827961179874,142.471471673523,200.990581616883,75.2439306555135,10.434727086962994,8.526633932683126,5.509135247188169,0.0692527749868838,-0.02355200000000046,24.315822490900985,59.15796117987399,-1.5385283264769782,47.540581616883,-17.416069344486502,-0.543272913037006,-0.4593660673168749,0.442135247188169,0.0133027749868838]


MgO_LewisCatlow= OrderedDict()
MgO_LewisCatlow['chrg_Mg'] = +2.0
MgO_LewisCatlow['chrg_O']  = -2.0
MgO_LewisCatlow['MgMg_A']   = 0.0 
MgO_LewisCatlow['MgMg_rho'] = 0.5
MgO_LewisCatlow['MgMg_C']   = 0.0
MgO_LewisCatlow['MgO_A']    = 821.6
MgO_LewisCatlow['MgO_rho']  = 0.3242
MgO_LewisCatlow['MgO_C']    = 0.0
MgO_LewisCatlow['OO_A']     = 2274.00 
MgO_LewisCatlow['OO_rho']   = 0.1490
MgO_LewisCatlow['OO_C']     = 27.88

MgO_BallGrimes1 = OrdredDict()
MgO_BallGrimes1['chrg_Mg'] = 2.0
MgO_BallGrimes1['chrg_O'] = -2.0
MgO_BallGrimes1['MgMg_A'] = 0.0
MgO_BallGrimes1['MgMg_rho'] 0.5
MgO_BallGrimes1['MgMg_C'] = 0.0
MgO_BallGrimes1['MgO_A'] = 1279.69
MgO_BallGrimes1['MgO_rho'] = 0.29969
MgO_BallGrimes1['MgO_C'] =0.0
MgO_BallGrimes1['OO_A'] = 9547.96
MgO_BallGrimes1['OO_rho'] = 0.21916
MgO_BallGrimes1['OO_C'] = 32.0

MgO_BallGrimes2 = OrderedDict()
MgO_BallGrimes2['chrg_Mg'] = 1.7
MgO_BallGrimes2['chrg_O'] = -1,7
MgO_BallGrimes2['MgMg_A'] = 0.0
MgO_BallGrimes2['MgMg_rho'] = 0.5
MgO_BallGrimes2['MgMg_C'] = 0.0
MgO_BallGrimes2['MgO_A'] = 929.69
MgO_BallGrimes2['MgO_rho'] 0.29909
MgO_BallGrimes2['MgO_C'] = 0.0
MgO_BallGrimes2['OO_A' = 4870
MgO_BallGrimes2['OO_rho'] = 0.2679
MgO_BallGrimes2['OO_C'] = 77
pyposmat_optimal_filename ='optimal_1000.out'

def create_pyposmat_optimal_file(
        filename_in,
        filename_out,
        n,
        scaling_factors,
        qoi_references,
        err_type):
    assert type(filename_in) == str
    assert type(filename_out) == str
    assert n == int
    assert scaling_factors == str
    assert isinstance(qoi_references,dict)
    assert scaling_factors in qoi_references
    assert err_type == str

    obj_datafile = PyposmatDataFile(filename=filename_in)
    obj_datafile.read(filename=filename_in)

    obj_datafile.qoi_references = copy.deepcopy(qoi_references)
    obj_datafile.create_optimal_population(
            n=n,
            scaling_factors=scaling_factors,
            err_type=err_type)

    obj_datafile.write_optimal_population(
            filename=filename_out,
            n=n,
            scaling_factors=scaling_factors,
            err_type=err_type)

if __name__ == "__main__":
    n = 1000
    scaling_factors='DFT'
    err_type='abs'
    pypospack_optimal_filename='optimal_{n}.out'.format(n=n)
    is_create_file = False

    if is_create_file:
        qoi_references = OrderedDict()
        qoi_references['DFT'] = copy.deepcopy(qoi_reference_dft)
        create_pypospack_optimal_file(
                filename_in=pypospack_data_filename,
                filename_out=pypospack_optimal_filename,
                n = 1000,
                scaling_factors='DFT',
                qoi_references = copy.deepcopy(qoi_references),
                err_type='abs')

    #--------------------------------------------------------------------------
    # PHONON EXAMPLE
    #--------------------------------------------------------------------------
    from pypospack.task.gulp import GulpPhononCalculation
    from pypospack.task.gulp import GulpGammaPointPhonons
    MgO_phonon_task= OrderedDict()
    MgO_phonon_task['task_name'] = 'MgO_NaCl.phonon'
    MgO_phonon_task['task_type'] = 'gulp_gamma_phonons'
    MgO_phonon_task['task_directory'] = 'MgO_NaCl.phonon'
    MgO_phonon_task['structure_filename'] = os.path.join(
            'test_PyposmatFileSampler',
            'MgO_NaCl_prim.gga.relax.vasp')
   
    MgO_potential = OrderedDict()
    MgO_potential['potential_type'] = 'buckingham'
    MgO_potential['symbols'] = ['Mg','O']
    MgO_potential['parameter_names'] = [
        'chrg_Mg','chrg_O',
        'MgMg_A','MgMg_rho','MgMg_C',
        'MgO_A','MgO_rho','MgO_C',
        'OO_A','OO_rho','OO_C']
    
    MgO_configuration_phonon = OrderedDict()
    MgO_configuration_phonon['potential'] \
            = copy.deepcopy(MgO_potential)
    MgO_configuration_phonon['tasks'] = OrderedDict()
    MgO_configuration_phonon['tasks']['MgO_NaCl.phonon'] \
            = copy.deepcopy(MgO_phonon_task)
    MgO_configuration_phonon['parameters'] = None
    #<--- use generic variables so we can test different stuff
    configuration = copy.deepcopy(MgO_configuration_phonon)
    task_information = copy.deepcopy(MgO_phonon_task)
    #<---- code start
    file_sampler = PyposmatFileSampler(
            filename_in = pypospack_optimal_filename)
    file_sampler.add_task(
            task_name=task_information['task_name'],
            task_directory=task_information['task_directory'],
            structure_filename=task_information['structure_filename'],
            task_type=task_information['task_type'])
    #<--- testing
    _task_name = task_information['task_name']
    assert _task_name in file_sampler.tasks
    assert type(file_sampler.tasks[_task_name])== GulpGammaPointPhonons
    assert file_sampler.tasks[_task_name].status == 'INIT'
    #<--- end testing
    file_sampler.configuration = copy.deepcopy(configuration)

    #<--- evaluate just one
    file_sampler.evaluate_parameter_set(MgO_LewisCatlow)
    for task_name,task in file_sampler.tasks.items(): 
        print(task.results)
    #<--- evaluate MgO
    #file_sampler.run()
    #file_sampler.read_pyposmat_datafile()
    #file_sampler.sample_from_file()
    #print('n_rows={}'.format(file_sampler.parameter_df.shape[0]))
