mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.7198127620048314
set group O charge -1.7198127620048314

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1134.8750919310144 0.3297433050145011 0.0 ${R_cut}
pair_coeff 2 2 21118.401484434533 0.1794796080547274 62.6057278667615 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
