mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.6546306328410598
set group O charge -1.6546306328410598

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 983.06062191212 0.30106964934880476 0.0 ${R_cut}
pair_coeff 2 2 1816.387023346269 0.2483974351081434 66.76648572299771 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
