import ase.build
import numpy as np


class DiamondStructure(object):

    def __init__(self):
        pass

    def build_from_ase(self):
        ase_cell = bulk('Si',5.431,cubic=True)

def print_ase_cell_details(ase_cell):
    print(ase_cell.cell)
    for atom in ase_cell:
        print(atom)

def make_silicon_bulk():

    ase_cell = ase.build.bulk('Si',a=5.431, cubic=True)
    return ase_cell

def make_silicon_supercell():
    print(80*'=')
    print('{:^80}'.format('making ASE supercell'))
    print(80*'=')

    symbols=['Si']
    a=5.431
    sc=[3,3,3]

    if isinstance(sc,list):
        P = np.matmul(np.identity(3),np.array(sc))
        print(P)
    elif isinstance(sc,np.ndarray):
        P = np.matmul(np.identity(3),sc)

    ase_cell = ase.build.make_supercell(
            prim=ase.build.bulk(''.join(symbols),a=a,cubic=True),
            P=np.identity(3)*np.array(sc)
    )
    return ase_cell

def get_interatomic_distance(ase_cell,i,j):
    _,r_i = ase_cell.get_positions[i]
    _,r_j = ase_cell.get_positions[j]
    
def find_closest_atom_to_coordinate(ase_cell,x):
    if isinstance(x,list):
        r_i = np.array(x)
    elif isintance(x,np.ndarray):
        r_i = x
    else:
        raise TypeError()

    all_r_j = [atom[1] for atom in ase_cell]
    

def make_monovacancy_silicon_supercell():
    print(80*'=')
    print('{:^80}'.format('making monovacancy structure'))
    print(80*'=')

    symbols=['Si']
    a=5.431
    sc=[3,3,3]

    if isinstance(sc,list):
        P = np.matmul(np.identity(3),np.array(sc))
        print(P)
    elif isinstance(sc,np.ndarray):
        P = np.matmul(np.identity(3),sc)

    find_closest_atom_to_coordinate(ase_cell,x)    
    ase_cell = ase.build.make_supercell(
            prim=ase.build.bulk(''.join(symbols),a=a,cubic=True),
            P=np.identity(3)*np.array(sc)
    )
    print_ase_cell_details(ase_cell=ase_cell)


if __name__ == "__main__":
    si_unit = make_silicon_bulk()
    #print_ase_cell_details(ase_cell=si_bulk)
    print(si_unit.cell)
    print(si_unit.get_positions())
    make_silicon_supercell()
