mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.637051747189104
set group O charge -1.637051747189104

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 855.6693654784183 0.31766411832708447 0.0 ${R_cut}
pair_coeff 2 2 11792.771131333073 0.1362741653681459 36.4169906289783 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
