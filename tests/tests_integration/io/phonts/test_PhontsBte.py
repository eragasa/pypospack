import pytest

def test__import__from_pypospack_task_phonts():
    from pypospack.task.phonts import PhontsBte


def test__import__from_pypospack_io_phonts():
    from pypospack.io.phonts import PhontsInputFile

def test_PhontsBte():
    from pypospack.task.phonts import PhontsBte

    phonts_bte = PhontsBte(
            task_name='phonts',
            task_directory='phonts',
            phonts_filename='phonon_input.dat',
            structure_filename='Si_dia_unit.relax.gga.vasp',
            restart=True)




