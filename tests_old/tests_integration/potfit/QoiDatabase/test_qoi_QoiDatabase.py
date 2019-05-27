import pypospack.potfit as potfit

class TestQoiDatabase(object):

    def setup__init_wo_yaml_file(self):
        self.QoiDatabase = potfit.QoiDatabase()

    def test__init_wo_yaml_file(self):
        self.setup__init_wo_yaml_file()
        assert self.QoiDatabase.filename == "pypospack.qoi.yaml"
        assert isinstance(self.QoiDatabase.qois, dict)

    def test_read(self):
        self.setup__init_wo_yaml_file()
        self.QoiDatabase.read("pypospack.qoi.yaml")

        assert self.QoiDatabase.filename == "pypospack.qoi.yaml"
        assert isinstance(self.QoiDatabase.qois, dict)
    
    def test_add_qoi(self):
        name = 'name'
        qoi = 'qoi'
        structures = ['structure']
        target = 1.0

        self.setup__init_wo_yaml_file()
        self.QoiDatabase.add_qoi(name,qoi,structures,target)

        assert self.QoiDatabase.qois[name]['qoi'] == qoi
        assert set(self.QoiDatabase.qois[name]['structures']) == set(structures)
        assert self.QoiDatabase.qois[name]['target'] == target
    #--------------------------------------------------------------    
    def setup__init_w_yaml_file(self):
        self.QoiDatabase = potfit.QoiDatabase(filename = "pypospack.qoi.yaml")

    def test_QoiDatabase__init_w_yaml_file(self):
        self.setup__init_w_yaml_file()
         
        assert self.QoiDatabase.filename == "pypospack.qoi.yaml"
        assert isinstance(self.QoiDatabase.qois, dict)
        

if __name__ == "__main__":
    qoi_info = potfit.QoiDatabase()
    qoi_info.read("pypospack.qoi.yaml")
    
    qoi_info.check()
    
    for k,v in qoi_info.qois.items():
        print(k,v)


