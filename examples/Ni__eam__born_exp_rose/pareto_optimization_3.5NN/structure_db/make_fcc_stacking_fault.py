import copy
import numpy as np
from pypospack.crystal import SimulationCell
from pypospack.io.ase import make_fcc_111_cell

def make_fcc_stacking_fault(symbols,stacking_sequence):
    # use ASE to make fcc 111 oriented cell
    ase_fcc_111 = make_fcc_111_cell(symbols=symbols)
    ase_lattice = ase_fcc_111.get_cell()

    for a in ase_fcc_111:
        for i in range(3):
            a.position[i]=round(a.position[i],8)

    z_layers = []
    for a in ase_fcc_111:
        if a.position[2] not in z_layers:
            z_layers.append(a.position[2])

    n_layers = len([v for v in stacking_sequence if v is not ','])
    layer_spacing = ase_lattice[2,2]/len(z_layers)

    simulation_cell=SimulationCell()
    simulation_cell.lattice_parameter = 1.0

    simulation_cell.H = copy.deepcopy(ase_lattice)
    simulation_cell.H[2,2] = layer_spacing*n_layers
    n_layers = len([v for v in stacking_sequence if v is not ','])
    for n,c in enumerate([v for v in stacking_sequence if v!=',']):
        if c == 'A':
            layer = 'A'
            z_layer = z_layers[0]
        elif c == 'B':
            layer = 'B'
            z_layer = z_layers[1]
        elif c == 'C':
            layer = 'C'
            z_layer = z_layers[2]
        for a in ase_fcc_111:
            if a.position[2] == z_layer:
                s  = a.symbol
                x  = copy.deepcopy(a.position)
                x[0] = x[0]/simulation_cell.H[0,0]
                x[1] = x[1]/simulation_cell.H[1,1]
                x[2] = n*layer_spacing/simulation_cell.H[2,2]
                simulation_cell.add_atom(s,x)
    return simulation_cell
if __name__ == "__main__":
    from pypospack.io.vasp import Poscar
    # INTRINSIC STACKING FAULT
    symbols = ['Ni']
    stacking_sequence = "ABC,ABC,AC,ABC,ABC"
    simulation_cell = make_fcc_stacking_fault(symbols,stacking_sequence)
    simulation_cell.normalize_h_matrix()
    poscar = Poscar(simulation_cell)
    poscar.write('Ni_fcc_isf.vasp')

    symbols = ['Ni']
    stacking_sequence = "ABC,ABC,ABAC,ABC,ABC"
    simulation_cell = make_fcc_stacking_fault(symbols,stacking_sequence)
    simulation_cell.normalize_h_matrix()
    poscar = Poscar(simulation_cell)
    poscar.write('Ni_fcc_esf.vasp')

    #symbols = ['Ni']
    #stacking_sequence = "ABC,ABC,AABC,ABC,ABC"
    #simulation_cell = make_fcc_stacking_fault(symbols,stacking_sequence)
    #simulation_cell.normalize_h_matrix()
    #poscar = Poscar(simulation_cell)
    #poscar.write('Ni_fcc_usf.vasp')
