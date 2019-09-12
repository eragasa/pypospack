import os
from pypospack.crystal import SimulationCell
from pypospack.crystal import FaceCenteredCubic
from pypospack.io.vasp.simulation import VaspSimulation
from pypospack.dft import determine_kpoint_meshes

slurm_default = [

]
def initialize_poscar(cell):
    if isinstance(cell, str):
        poscar = Poscar()
        poscar.read(filename=cell)
    elif isinstance(cell, SimulationCell):
        poscar = Poscar(obj_cell=cell)

    return poscar

def get_base_simulation():
    base_simulation = VaspSimulation()
    base_simulation.poscar = initialize_poscar(cell)
    return base_simulation

kpoint_meshes = determine_kpoint_meshes(simulation_cell=fcc_cell)
for k,v in kpoint_meshes.items():
    simulation_path = os.path.join(simulation_base_path,k)
    print(simulation_path)
