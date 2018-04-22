mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.32779280783179
set group O charge -2.32779280783179

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 996.6612573293085 0.31791195727041904 0.0 ${R_cut}
pair_coeff 2 2 3586.835772764404 0.16274043818444495 33.475942371618956 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
