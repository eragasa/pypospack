mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.66768155913
set group O charge -1.66768155913

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1248.0823356 0.293822485762 0.0 ${R_cut}
pair_coeff 2 2 31778.9957117 0.212502502632 34.6036277069 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
