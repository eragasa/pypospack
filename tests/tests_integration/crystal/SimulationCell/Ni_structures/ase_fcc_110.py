import ase.build

atoms = ase.build.fcc110('Ni',size=(1,1,1),pbc=(True,True,True))
print(atoms.cell)
print(atoms)
