mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.1834714935002904
set group O charge -2.1834714935002904

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 950.3965262142411 0.29958011010557994 0.0 ${R_cut}
pair_coeff 2 2 14050.592204812488 0.2773879721194452 73.85408645436596 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
