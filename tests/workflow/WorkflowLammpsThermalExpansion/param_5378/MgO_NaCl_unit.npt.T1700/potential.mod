mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.67557694923
set group O charge -1.67557694923

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1405.50211521 0.282257384366 0.0 ${R_cut}
pair_coeff 2 2 19860.1058829 0.143671361936 6.84969713079 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
