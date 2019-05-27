import pytest
import os,shutil
from collections import OrderedDict

def test__import__from_pypospack_task_phonts():
    from pypospack.task.phonts import PhontsBte


def test__import__from_pypospack_io_phonts():
    from pypospack.io.phonts import PhontsInputFile

@pytest.fixture
def resource_phontsbte(request):
    def make_test_package(
            task_name,
            task_directory,
            phonts_filename,
            structure_filename,
            restart
    ):
        from pypospack.task.phonts import PhontsBte
        phonts_bte = PhontsBte(
                task_name=task_name,
                task_directory=task_directory,
                phonts_filename=phonts_filename,
                structure_filename=structure_filename,
                restart=restart)
        
        #<--- register the teardown
        def fin():
            print('....deleting the directory:{}'.format(task_directory))
            shutil.rmtree(task_directory)
        request.addfinalizer(fin)
        
        return phonts_bte

    return make_test_package

test_one = OrderedDict()
test_one['class_args']=OrderedDict()
test_one['class_args']['task_name']='phonts_name'
test_one['class_args']['task_directory']='phonts_directory'
test_one['class_args']['phonts_filename']='phonon_input.dat'
test_one['class_args']['structure_filename']=os.path.join('test_PhontsBte','Si_dia_unit.relax.gga.vasp')
test_one['class_args']['restart']=True


@pytest.mark.parametrize("class_args",[
    test_one['class_args']
])
def test_PhontsBte(class_args,resource_phontsbte):
    phonts_bte = resource_phontsbte(
            task_name=class_args['task_name'],
            task_directory=class_args['task_directory'],
            phonts_filename=class_args['phonts_filename'],
            structure_filename=class_args['structure_filename'],
            restart=class_args['restart'])

    #from pypospack.task.phonts import PhontsBte
    #phonts_bte = PhontsBte(
    #        task_name='phonts',
    #        task_directory='phonts',
    #        phonts_filename='phonon_input.dat',
    #        structure_filename='Si_dia_unit.relax.gga.vasp',
    #        restart=True)

if False:
    def test_PhontsBte():
        from pypospack.task.phonts import PhontsBte
        phonts_bte = PhontsBte(
                task_name='phonts',
                task_directory='phonts',
                phonts_filename='phonon_input.dat',
                structure_filename='Si_dia_unit.relax.gga.vasp',
                restart=True)




