import pypospack.potfit as potfit

def test__import__from_pypospack_qoi():
    from pypospack.qoi import CrystalStructureGeometry

def test____import__():
    from pypospack.qoi import CrystalStructureGeometry

    obj_qoi  = CrystalStructureGeometry(
            qoi_name = qoi_name,
            structures = qoi_structures)

    assert obj_qoi.qoi_name == qoi_name
    assert obj_qoi.qoi_type == qoi_type
    assert obj_qoi.required_simulations is None
    assert type(obj_qoi.required_tasks) == OrderedDict
    assert type(obj_qoi.required_variables) == OrderedDict
    assert type(obj_qoi.predicted_values) is None

if __name__ == "__main__":
    from pypospack.qoi import CrystalStructureGeometry
    qoi_structures = ['NaCl_unit']
    qoi_name = 'NaCl.geometry'

    obj_qoi = CrystalStructureGeometry(
            qoi_name = qoi_name,
            structures = qoi_structures)

    print(obj_qoi.get_required_simulations())
