import pypospack.potfit as potfit
import pypospack.qoi as qoi
class TestQoiManager(object):
    def test_without_QoiInfo(self):
        self.qoi_manager = QoiManager()
        assert self.qoi_manager.qoi_names is None

    def test_with_QoiInfo(self):
        self.qoi_info = potfit.QoiDatabase()
        self.qoi_info.read("pypospack.qoi.yaml")
        self.qoi_manager = QoiManager(self.qoi_info)
        assert isinstance(self.qoi_manager.qoi_names, list)
        assert isinstance(self.qoi_manager.supported_qoi, list)

if __name__ == "__main__":
    qoi_info = potfit.QoiDatabase()
    qoi_info.read("pypospack.qoi.yaml")
    qoi_manager = qoi.QoiManager(qoi_info)
   
    s = ( " This script determines what simulations what a yaml configured file "
          " will produce in pypospack potential optimization software ")
    print(s)
    print(80*'-')
    print('qoi_names')
    for q in qoi_manager.qoi_names:
        print('\t{}'.format(q))


    required_simulations = qoi_manager.get_required_simulations()
    print(80*'-')
    print('required_simulations:')
    for k,v in required_simulations.items():
        print(k)
        for kv, vv in v.items():
            print("\t{}:{}".format(kv,vv))
    
    # Here we have the problem of how to deal with prequisite simulations
    # being complete before the other one completes.  Here we borrow the
    # terminology of GANNT chart.

    # Task: unique_task_name -> structure_name.simulation
    # Predecessor: unique_task_name -> structure_name.simulation
    # Predecessor_Info:
    #     global_variable_name: structure_name.simulation_type.variable_name
    #     local_variable_name: variable_name
    #     variable_value: 
    print(80*'-')
    print('qois_objects:')
    for q,obj in qoi_manager.obj_qois.items():
        print('\t',q,type(obj))
