mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.5471876556322004
set group O charge -1.5471876556322004

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1187.8352669893875 0.2999606674467126 0.0 ${R_cut}
pair_coeff 2 2 19707.723693181426 0.36634245735287185 45.119127443553396 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
