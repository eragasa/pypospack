from pypospack.task import Task
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp

def poscar_to_gulp_string(poscar):
    if isinstance(poscar,crystal.SimulationCell):
        return simulation_cell_to_gulp_string(poscar)
    elif isinstance(poscar,str):
        simcell = vasp.Poscar()
        simcell.read(poscar)
        return simulation_cell_to_gulp_string(simcell)

def simulation_cell_to_gulp_string(sim_cell):
    if not isinstance(sim_cell,crystal.SimulationCell):
        raise ValueError('simcell is not an instance of pypospack simulation cell')

    H = sim_cell.H * sim_cell.a0
    str_out = "vectors\n"
    str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[0,0],H[0,1],H[0,2])
    str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[1,0],H[1,1],H[1,2])
    str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[2,0],H[2,1],H[2,2])
    str_out += "fractional\n"
    for s in sim_cell.symbols:
        for a in sim_cell.atomic_basis:
            if a.symbol == s:
                try:
                    str_out += "{} core {} {} {}\n".format(\
                            s,a.position[0],a.position[1],a.position[2])
                except:
                    print(s)
                    print(a.symbol)
                    print(a.position)
                    raise
    return str_out

class GulpSimulation(Task):
    def __init__(self,task_name,task_directory,restart=False):
        Task.__init__(self,task_name,task_directory,restart)
        self.structure = vasp.Poscar()
        self.potential = None
        self.status == 'INIT'

    def config(self,config_dict,param_dict):

   
