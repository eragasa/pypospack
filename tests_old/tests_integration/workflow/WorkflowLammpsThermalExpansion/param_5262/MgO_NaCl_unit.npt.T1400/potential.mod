mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.695400086
set group O charge -1.695400086

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1303.37543516 0.290482854228 0.0 ${R_cut}
pair_coeff 2 2 9936.88636401 0.22003654812 26.9596347598 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
