import pytest


def test__import__from_pypospack_qoi():
    from pypospack.qoi import Qoi

qoi_name='qoi_name'
qoi_type='qoi_type'
structure_names=['structure_filename']

def test____init__():
    from pypospack.qoi import Qoi

    qoi = Qoi(
            qoi_name=qoi_name,
            qoi_type=qoi_type,
            structure_names=structure_names)

    assert qoi.qoi_name == qoi_name
    assert qoi.qoi_type == qoi_type
    assert qoi.structure_names == structure_names
    assert qoi.reference_values is None
    assert qoi.task_definitions is None
    assert qoi.task_dependencies is None
    assert qoi.results is  None
