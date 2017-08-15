import pypospack.io.vasp as vasp

poscar_filename_in = "Ni_fcc_unit.vasp"
poscar = vasp.Poscar()
poscar.read(filename=poscar_filename_in)

print(poscar.comment)
print(poscar.a1)
print(poscar.a2)
print(poscar.a3)
print(poscar.h_matrix[0,:])
print(poscar.h_matrix[1,:])
print(poscar.h_matrix[2,:])
print(len(poscar.atoms))
for a in poscar.atoms:
    print(a.symbol,a.position,a.magnetic_moment)
