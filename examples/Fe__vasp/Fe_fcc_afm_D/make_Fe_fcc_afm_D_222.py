import os, shutil, subprocess, copy
from pypospack.crystal import SimulationCell
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as minimize_init_slurm

def make_bulk_structure(filename_poscar_in,
                        filename_outcar_in=None,
                        supercell=None):
    poscar_in = vasp.Poscar()
    poscar_in.read(filename=filename_poscar_in)

    if filename_outcar_in is not None:
        # add magnetic moments to the structure
        outcar_in = vasp.Outcar()
        outcar_in.read(filename=filename_outcar_in)
        magnetic_moments = outcar_in.get_magnetic_moments()
        for i, atom in enumerate(poscar_in.atomic_basis):
            atom.magmom = magnetic_moments[i][-1]
        z_distances = [k.position[2] for k in poscar_in.atomic_basis]
        min_z_distance = min(z_distances)
        for atom in poscar_in.atomic_basis:
            atom.position[2] = atom.position[2] - min_z_distance

    # make supercell
    if supercell is not None:
        poscar_out = vasp.Poscar(
            crystal.make_super_cell(
                structure=poscar_in,
                sc=supercell
            )
        )
    else:
        poscar_out = poscar_in

    # sort by z_axis
    new_atomic_basis = []
    while len(poscar_out.atomic_basis) > 0:
        z_distances = [k.position[2] for k in poscar_out.atomic_basis]
        idx_min = z_distances.index(min(z_distances))
        new_atomic_basis.append(copy.deepcopy(poscar_out.atomic_basis[idx_min]))
        del poscar_out.atomic_basis[idx_min]
    poscar_out.atomic_basis = new_atomic_basis
    return poscar_out
    #z_distances = [k.position[2] for k in poscar_out.atomic_basis]
    #poscar_out.write(system_name+".init.vasp")

minimization_path = "minimization"
system_name = "Fe_fcc_afm_D_222"
simulation_directory = system_name
filename_poscar_in = os.path.join(minimization_path,'CONTCAR')
filename_outcar_in = os.path.join(minimization_path,'OUTCAR')
supercell = [4,4,2]

structures = {
    'Fe_fcc_afm_D':[],
    'Fe_fcc_afm_D_vac_A':
        [{
            'type':'vac',
            'symbol':'Fe',
            'position':[0.250,0.000,0.50000]
         }],
    'Fe_fcc_afm_D_vac_B':
        [{
            'type':'vac',
            'symbol':'Fe',
            'position':[0.125,0.000,0.61976]
         }],
    'Fe_fcc_afm_D_vac_C':
        [{
            'type':'vac',
            'symbol':'Fe',
            'position':[0.125,0.125,0.50000]
         }],
    'Fe_fcc_afm_D_vac_D':
        [{
            'type':'vac',
            'symbol':'Fe',
            'position':[0.125,0.000,0.36976]
         }],
    'Fe_fcc_afm_D_vac_E':
        [{
            'type':'vac',
            'symbol':'Fe',
            'position':[0.000,0.000,0.75000]
         }],
    'Fe_fcc_afm_D_vac_F':
        [{
            'type':'vac',
            'symbol':'Fe',
            'position':[0.250,0.125,0.61965]
         }]
}

for name,defects in structures.items():
    poscar_out = make_bulk_structure(filename_poscar_in=filename_poscar_in,
                                     filename_outcar_in=filename_outcar_in,
                                     supercell=supercell)
    filename_poscar_out = "{}.init.vasp".format(name)
    print(filename_poscar_out)

    for defect in defects:
        defect_type = defect['type']
        if defect_type == 'vac':
            defect_symbol = defect['symbol']
            defect_position = defect['position']
            poscar_out.add_vacancy(symbol=defect_symbol,
                                   position=defect_position)
        elif defect_type == 'interstitial':
            defect_symbol = defect['symbol']
            defect_position = defect['position']
            try:
                defect_magmom = defect['magmom']
            except KeyError as e:
                defect_magmom = 0.0
            poscar_out.add_interstitial(symbol,position,magmom)
    poscar_out.write(filename_poscar_out)

    with open("{}.init.magmom".format(name),'w') as f:
        print('\t'+poscar_out.get_magmom_tag())
        f.write(poscar_out.get_magmom_tag())
