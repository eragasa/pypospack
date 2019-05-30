import pytest

from collections import OrderedDict

def test__init__with_None():
    _qoidb = None
    from pypospack.qoi import QoiDatabase
    from pypospack.qoi import QoiManager
    qoimanager = QoiManager(qoi_database=_qoidb,fullauto=False)

    assert isinstance(qoimanager.qoidb,QoiDatabase)
    #assert isinstance(qoimanager.qoidb.qoi_names,list)
    #assert isinstance(qoimanager,qoi_names,list)

def test__init__with_implicit_None():
    from pypospack.qoi import QoiDatabase
    from pypospack.qoi import QoiManager
    qoimanager = QoiManager(fullauto=False)

    assert isinstance(qoimanager.qoidb,QoiDatabase)
    #assert isinstance(qoimanager.qoidb.qoi_names,list)
    #assert isinstance(qoimanager,qoi_names,list)

if True:
    def test__init__with_QoiDatabase():
        _qoidb_filename_in = "pypospack.qoi.yaml"

        from pypospack.qoi import QoiDatabase
        qoidb = QoiDatabase()
        qoidb.read(filename=_qoidb_filename_in)

        from pypospack.qoi import QoiManager
        qoimanager = QoiManager(qoi_database=qoidb,fullauto=False)

        assert isinstance(qoimanager.qoidb,QoiDatabase)
        assert isinstance(qoimanager.qoidb.qoi_names,list)
        assert isinstance(qoimanager.qoi_names,list)

    def test__init__with_filename():
        _qoidb_filename_in = "pypospack.qoi.yaml"

        from pypospack.qoi import QoiDatabase
        from pypospack.qoi import QoiManager
        qoimanager = QoiManager(qoi_database=_qoidb_filename_in,fullauto=False)

        assert isinstance(qoimanager.qoidb,QoiDatabase)
        assert isinstance(qoimanager.qoidb.qoi_names,list)
        assert isinstance(qoimanager.qoi_names,list)

    def test__configure():
        _qoidb_filename_in = "pypospack.qoi.yaml"

        from pypospack.qoi import QoiDatabase
        from pypospack.qoi import QoiManager
        qoimanager = QoiManager(qoi_database=_qoidb_filename_in,fullauto=False)
        qoimanager.configure()

        assert isinstance(qoimanager.obj_Qoi,OrderedDict)
        assert all([isinstance(v,qoi.Qoi) for k,v in qoimanager.obj_Qoi.items()])

    def test__get_required_tasks():
        _qoidb_filename_in = "pypospack.qoi.yaml"

        from pypospack.qoi import QoiDatabase
        from pypospack.qoi import QoiManager
        qoimanager = QoiManager(qoi_database=_qoidb_filename_in,fullauto=False)
        qoimanager.configure()
        qoimanager.determine_required_tasks()
        
        assert isinstance(qoimanager.required_tasks,OrderedDict)


