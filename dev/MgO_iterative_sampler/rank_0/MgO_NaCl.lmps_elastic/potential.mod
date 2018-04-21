mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.467095066210368
set group O charge -2.467095066210368

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 935.4807786808233 0.3084619384680258 0.0 ${R_cut}
pair_coeff 2 2 12266.368986574063 0.11126335340408924 33.77611037763418 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
thermo	1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
