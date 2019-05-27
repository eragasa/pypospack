# ---- initialize simulations
clear
units metal
dimension 3
boundary p p p
atom_style atomic
atom_modify map array
# ---- create atoms
read_data lammps.structure
# ---- define interatomic potential
pair_style eam/alloy
pair_coeff * * eam.alloy Ni
neighbor 2.0 bin
neigh_modify delay 10 check yes

# ---- define settings
compute eng all pe/atom
compute eatoms all reduce sum c_eng

# ---- run minimization
reset_timestep 0
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms

run 0

variable natoms equal "count(all)"
variable tot_energy equal "c_eatoms"
variable length_x equal "lx/4"
variable length_y equal "ly/4"
variable length_z equal "lz/4"
variable ecoh equal "pe/atoms"
variable press_t equal "press"
variable press_x equal "pxx"
variable press_y equal "pyy"
variable press_z equal "pzz"
#variable ecoh equal "v_etotal/v_atoms"
# --- output
print "pyPosMat output section"
print "tot_energy = ${tot_energy}"
print "num_atoms = ${natoms}"
print "latt_const_a = ${length_x}"
print "latt_const_b = ${length_y}"
print "latt_const_c = ${length_z}"
print "pressure_total = ${press_t}"
print "pressure_xx = ${press_x}"
print "pressure_yy = ${press_y}"
print "pressure_zz = ${press_z}"
print "ecoh = ${ecoh}"
print "lammps_sim_done"
