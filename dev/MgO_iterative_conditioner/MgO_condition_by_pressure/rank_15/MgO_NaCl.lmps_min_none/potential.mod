mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.8710163124128623
set group O charge -1.8710163124128623

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1144.968548749854 0.300234270752247 0.0 ${R_cut}
pair_coeff 2 2 17151.24735745269 0.3396416962302238 53.17683019667447 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
