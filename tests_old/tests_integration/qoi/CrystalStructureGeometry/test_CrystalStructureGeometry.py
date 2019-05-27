import pytest
from collections import OrderedDict

MgO_expected_results = OrderedDict()
MgO_expected_results['qoi_type'] = 'lmps_relax_all'
MgO_expected_results['structure_ideal'] = 'MgO_NaCl'

MgO_qoi_configuration1 = OrderedDict()
MgO_qoi_configuration1['qoi_name'] = 'MgO_NaCl.lmps_relax_all'
MgO_qoi_configuration1['structures'] = 'MgO_NaCl'
MgO_qoi_configuration1['expected'] = MgO_expected_results

MgO_qoi_configuration2 = OrderedDict()
MgO_qoi_configuration2['qoi_name'] = 'MgO_NaCl.lmps_relax_all'
MgO_qoi_configuration2['structures'] = OrderedDict()
MgO_qoi_configuration2['structures']['ideal'] = 'MgO_NaCl'
MgO_qoi_configuration2['expected'] = MgO_expected_results

MgO_qoi_configuration3 = OrderedDict()
MgO_qoi_configuration3['qoi_name'] = 'MgO_NaCl.lmps_relax_all'
MgO_qoi_configuration3['structures'] = OrderedDict()
MgO_qoi_configuration3['structures'] = ['MgO_NaCl']
MgO_qoi_configuration3['expected'] = MgO_expected_results

testdata_id = [
        'type(MgO_structure) == str',
        'type(MgO_structure) == OrderedDict',
        'type(MgO_structure) == list'
    ]
testdata = [
        MgO_qoi_configuration1,
        MgO_qoi_configuration2,
        MgO_qoi_configuration3,
    ]


def test__import__from_pypospack_qoi():
    from pypospack.qoi import RelaxedStructureCalculations

@pytest.mark.parametrize("configuration",testdata,ids=testdata_id)
def test____init__(configuration):
    from pypospack.qoi import RelaxedStructureCalculations
   
    _qois_calculated = ['a0_min_all','a1_min_all','a2_min_all','a3_min_all',
            'Ecoh_min_all']
    _qoi_name = configuration['qoi_name']
    _qoi_type = configuration['expected']['qoi_type']
    _structures = configuration['structures']
    _structure_ideal = configuration['expected']['structure_ideal']
    qoi = RelaxedStructureCalculations(
            qoi_name=_qoi_name,
            structures=_structures)
    assert all([q in qoi.qois_calculated for q in _qois_calculated])
    assert qoi.qoi_name == _qoi_name
    assert qoi.qoi_type == 'lmps_relax_all'
    assert isinstance(qoi.structures,OrderedDict)
    assert qoi.structures['ideal'] == _structure_ideal

def test__determine_required_simulations():
    _qoi_name = MgO_qoi_configuration1['qoi_name']
    _qoi_type = MgO_qoi_configuration1['expected']['qoi_type']
    _structures = MgO_qoi_configuration1['structures']
    _structure_ideal = MgO_qoi_configuration1['expected']['structure_ideal']
    
    from pypospack.qoi import RelaxedStructureCalculations
    qoi = RelaxedStructureCalculations(
            qoi_name=_qoi_name,
            structures=_structures)

    qoi.determine_required_simulations()
    
    _ideal_task_type = 'lmps_min_all'
    _ideal_task_name = '{}.lmps_min_all'.format(_structure_ideal)
    _ideal_task_structure = _structure_ideal
    assert isinstance(qoi.tasks,OrderedDict)
    print(qoi.tasks)
    assert _ideal_task_name in qoi.tasks
    assert qoi.tasks[_ideal_task_name]['task_type'] == _ideal_task_type
    assert qoi.tasks[_ideal_task_name]['task_structure'] == _ideal_task_structure

if __name__ == '__main__':
    pass
