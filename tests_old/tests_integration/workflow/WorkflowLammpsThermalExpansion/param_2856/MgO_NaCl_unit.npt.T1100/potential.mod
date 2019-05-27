mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.71287520642
set group O charge -1.71287520642

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1268.25800143 0.292150986054 0.0 ${R_cut}
pair_coeff 2 2 34972.2762623 0.190168709461 15.9037418312 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
