mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.217493753944798
set group O charge -2.217493753944798

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1298.7491517563526 0.31679719255367117 0.0 ${R_cut}
pair_coeff 2 2 15659.281989980458 0.35765115098044997 52.21000040425842 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
