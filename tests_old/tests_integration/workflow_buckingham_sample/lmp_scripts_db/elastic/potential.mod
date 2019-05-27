pair_style eam/alloy
pair_coeff * * eam.alloy Ni
neighbor 2.0 bin
neigh_modify delay 10 check yes

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
