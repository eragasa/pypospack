mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.66679564584
set group O charge -1.66679564584

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1241.53975785 0.290586946924 0.0 ${R_cut}
pair_coeff 2 2 7336.97641517 0.111647641293 25.3936876208 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
