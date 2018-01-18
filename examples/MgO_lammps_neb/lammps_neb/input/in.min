# ---- initialize simulations
clear
units metal
dimension 3
boundary p p p
atom_style charge
atom_modify map array

# ---- BULK MINIMIZATION
# ---- create atoms
read_data MgO_NaCl_333.structure
# ---- define interatomic potential
include potential.mod
# include modify_atom_info.mod

# ---- define settings
compute eng all pe/atom
compute eatoms all reduce sum c_eng

# ---- run minimization
reset_timestep 0
fix 1 all box/relax iso 0.0 vmax 0.001
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
min_style cg
minimize 1e-25 1e-25 5000 10000

variable bulk_natoms equal "count(all)"
variable bulk_tot_energy equal "c_eatoms"
variable bulk_length_x equal "lx"
variable bulk_length_y equal "ly"
variable bulk_length_z equal "lz"
variable bulk_ecoh equal "pe/atoms"
# --- output
print "pyPosMat output section"
print "tot_energy = ${bulk_tot_energy}"
print "num_atoms = ${bulk_natoms}"
print "latt_const_a = ${bulk_length_x}"
print "latt_const_b = ${bulk_length_y}"
print "latt_const_c = ${bulk_length_z}"
print "ecoh = ${bulk_ecoh}"

#------------------------------------------------------------------------------
# INIT IMAGE MINIMIZATION
#------------------------------------------------------------------------------

# ---- initialize simulations
clear
units metal
dimension 3
boundary p p p
atom_style charge
atom_modify map array

# --- create atoms
read_data MgO_NaCl_333_fr_a_0.structure
change_box all x final 0 ${bulk_length_x} y final 0 ${bulk_length_y} z final 0 ${bulk_length_z} remap  
# ---- define interatomic potential
include potential.mod
# include modify_atom_info.mod

# ---- define settings
compute eng all pe/atom
compute eatoms all reduce sum c_eng

# ---- run minimization
reset_timestep 0
#fix 1 all box/relax iso 0.0 vmax 0.001
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
min_style cg
minimize 1e-25 1e-25 5000 10000

variable init_natoms equal "count(all)"
variable init_tot_energy equal "c_eatoms"
variable init_length_x equal "lx"
variable init_length_y equal "ly"
variable init_length_z equal "lz"
variable init_ecoh equal "pe/atoms"
# --- output
print "pyPosMat output section"
print "tot_energy = ${init_tot_energy}"
print "num_atoms = ${init_natoms}"
print "latt_const_a = ${init_length_x}"
print "latt_const_b = ${init_length_y}"
print "latt_const_c = ${init_length_z}"
print "ecoh = ${init_ecoh}"
print "lammps_sim_done"

# --- write structure out
write_data init.MgO_fr_a nocoeff
#------------------------------------------------------------------------------
# --- FINAL IMAGE MINIMIZATION
#------------------------------------------------------------------------------
# ---- initialize simulations
clear
units metal
dimension 3
boundary p p p
atom_style charge
atom_modify map array

# --- create atoms
read_data MgO_NaCl_333_fr_a_1.structure
change_box all x final 0 ${bulk_length_x} y final 0 ${bulk_length_y} z final 0 ${bulk_length_z} remap  

# ---- define interatomic potential
include potential.mod
# include modify_atom_info.mod
# ---- define settings
compute eng all pe/atom
compute eatoms all reduce sum c_eng

# ---- run minimization
reset_timestep 0
#fix 1 all box/relax iso 0.0 vmax 0.001
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
min_style cg
minimize 1e-25 1e-25 5000 10000

variable natoms equal "count(all)"
variable tot_energy equal "c_eatoms"
variable length_x equal "lx"
variable length_y equal "ly"
variable length_z equal "lz"
variable ecoh equal "pe/atoms"
#variable ecoh equal "v_etotal/v_atoms"
# --- output
print "pyPosMat output section"
print "tot_energy = ${tot_energy}"
print "num_atoms = ${natoms}"
print "latt_const_a = ${length_x}"
print "latt_const_b = ${length_y}"
print "latt_const_c = ${length_z}"
print "ecoh = ${ecoh}"

# --- write structure_out
write_data final.MgO_fr_a nocoeff
