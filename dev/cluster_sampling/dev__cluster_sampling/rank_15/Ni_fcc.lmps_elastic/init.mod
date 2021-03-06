# ---- init.mod file
variable up equal 1.0e-6
units metal
dimension 3
boundary p p p
atom_style atomic
atom_modify map array
variable cfac equal 1.0e-4
variable cunits string GPa
# ---- define minimization parameters
variable etol equal 0.0
variable ftol equal 1.0e-10
variable maxiter equal 100
variable maxeval equal 1000
variable dmax equal 1.0e-2
# --- read data structure
read_data lammps.structure
