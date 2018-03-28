mass 1 28.086

group Si type 1

pair_style sw
pair_coeff * * Si.parameters Si


neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
thermo	1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
