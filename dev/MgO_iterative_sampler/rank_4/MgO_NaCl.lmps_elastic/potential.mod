mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.1466418992637664
set group O charge -2.1466418992637664

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1089.4134673757053 0.31058320000297174 0.0 ${R_cut}
pair_coeff 2 2 18983.219085212437 0.12673144022038657 73.4670102644071 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
thermo	1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
