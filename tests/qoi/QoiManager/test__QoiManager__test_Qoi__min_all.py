import pytest
from collections import OrderedDict
#import pypospack.potfit as potfit
import pypospack.qoi as qoi

qois = ['Ecoh_min_all','a1_min_all','a2_min_all','a3_min_all',
        'a11_min_all','a12_min_all','a13_min_all',
        'a21_min_all','a22_min_all','a23_min_all',
        'a31_min_all','a32_min_all','a33_min_all']

_qoidb_OrderedDict = OrderedDict()
_qoidb_OrderedDict['qoi_name'] = 'MgO_NaCl.Ecoh'
_qoidb_OrderedDict['qoi_type'] = 'Ecoh_min_all'
_qoidb_OrderedDict['structures'] = OrderedDict()
_qoidb_OrderedDict['ideal'] = 'MgO_NaCl'

@pytest.mark.parametrize("qoi",qois)
def test__get_required_tasks(qoi):
    _qoidb_filename_in = "pypospack.qoi.yaml"

    from pypospack.qoi import QoiDatabase
    from pypospack.qoi import QoiManager

    _qoi_type = qoi
    _qoi_structure = 'MgO_NaCl'
    _qoi_name = "{}.{}".format(_qoi_structure,_qoi_type)
    _qoi_structures = OrderedDict()
    _qoi_structures['ideal'] = 'MgO_NaCl'
    _qoidb_QoiDatabase = QoiDatabase()
    _qoidb_QoiDatabase.add_qoi(
            qoi_name = _qoi_name,
            qoi_type = _qoi_type,
            structures = _qoi_structures,
            target = 4.5)


    qoimanager = QoiManager(qoi_database=_qoidb_QoiDatabase,fullauto=False)
    qoimanager.configure()
    qoimanager.determine_tasks()

    _task_type = 'lmps_min_all'
    _qoi_task_name = "{}.{}".format(_qoi_structures['ideal'],_task_type)
    assert isinstance(qoimanager.tasks,OrderedDict)
    assert _qoi_task_name in qoimanager.tasks
    assert qoimanager.tasks[_qoi_task_name]['task_type'] == _task_type
    assert qoimanager.tasks[_qoi_task_name]['task_structure'] \
            == _qoi_structures['ideal']

