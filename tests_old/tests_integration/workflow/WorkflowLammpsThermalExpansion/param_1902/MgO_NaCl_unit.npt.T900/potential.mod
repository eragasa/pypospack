mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.69754800024
set group O charge -1.69754800024

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1340.93810979 0.285979914483 0.0 ${R_cut}
pair_coeff 2 2 11096.6994948 0.145827000496 3.24328494634 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
