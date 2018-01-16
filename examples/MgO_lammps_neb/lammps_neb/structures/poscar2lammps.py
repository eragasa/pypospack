import pyflamestk.base as base
import pyflamestk.vasp as vasp
import pyflamestk.lammps as lammps

# this script converts a poscar file into a lammps structure file

fname_poscar = ['MgO_NaCl_333.vasp', 
                'MgO_NaCl_333_fr_a_0.vasp', 
                'MgO_NaCl_333_fr_a_1.vasp']
fname_lammps = ['MgO_NaCl_333.structure',
                'MgO_NaCl_333_fr_a_0.structure',
                'MgO_NaCl_333_fr_a_1.structure']
sym_order = ['Mg','O']

n_structures = len(fname_poscar)
for i in range(n_structures):
    poscar = vasp.Poscar()
    poscar.read_file(fname_in=fname_poscar[i])
    poscar.normalize_h_matrix()

    lmp_structure_file = lammps.StructureFile(obj=poscar)
    lmp_structure_file.write_file(fname_out=fname_lammps[i])

