import pypospack.io.vasp as vasp
import pypospack.io.lammps as lammps

if __name__ == '__main__':
    filename_in = 'rsrc/MgO_NaCl_unit.vasp'
    filename_out = 'structure.dat'
    vasp_structure = vasp.Poscar()
    vasp_structure.read(filename_in)
    lammps_structure = lammps.LammpsStructure(\
            obj=vasp_structure)
    print("a0",vasp_structure.a0,lammps_structure.a0)
    print("H",vasp_structure.H,lammps_structure.H)
    lammps_structure.write(filename=filename_out)
 


