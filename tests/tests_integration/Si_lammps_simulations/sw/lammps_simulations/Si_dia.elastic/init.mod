variable up equal 1.0e-6
 
# metal units, elastic constants in GPa
units		metal
dimension        3
boundary	p p p
atom_style charge
atom_modify map array
variable cfac equal 1.0e-4
variable cunits string GPa

# Define minimization parameters
variable etol equal 0.0 
variable ftol equal 1.0e-10
variable maxiter equal 100
variable maxeval equal 1000
variable dmax equal 1.0e-2

read_data lammps.structure
