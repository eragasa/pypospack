from pypospack.crystal import SimulationCell, make_super_cell
from pypospack.io.vasp import Poscar


import numpy as np
def brute_minimum_image_convention(r,H,a0):

    # Basic idea = triple loop over all nearby image, pick shortest relative vector.
    # This is horribly slow, and for very skewed cells, n has to become very large.
    result = None
    n = 2
    for i0 in range(-n, n+1):
        for i1 in range(-n, n+1):
            for i2 in range(-n, n+1):
                rp = r+np.dot(a0*H.T, [i0,i1,i2])
                d = np.linalg.norm(rp)
                if (result is None) or (result[1] > d):
                    result = (rp, d)
    return result[0]

poscar_filename = 'Ni_fcc_unit.gga.relaxed.vasp'
cell_poscar = Poscar()
cell_poscar.read(poscar_filename)
cell_poscar.normalize_h_matrix()
cell_sc = make_super_cell(structure=cell_poscar,sc=[3,3,3])
_positions = [a.position for a in cell_sc.atomic_basis]
_n_atoms = len(_positions)
min_distances = []
for i in range(1,_n_atoms):
    r1 = np.array(_positions[0])
    r2 = np.array(_positions[i])
    r = r2-r1
    r_min = brute_minimum_image_convention(r=r,H=cell_sc.H,a0=cell_sc.a0)
    d = np.linalg.norm(r_min)
    min_distances.append(d*cell_poscar.a0*cell_poscar.a1)
min_distances = np.array(min_distances)
distances = np.linspace(
        np.min(min_distances),
        np.max(min_distances),
        100)
np_lt_d = np.zeros(shape=distances.shape)
for i,d in enumerate(distances):
    print(d,np.where(min_distances <=d)[0].size)
    np_lt_d[i] = np.where(min_distances <= d)[0].size

print(distances)
print(np_lt_d)
