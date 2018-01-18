mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.73499569705
set group O charge -1.73499569705

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 969.401486338 0.314807375333 0.0 ${R_cut}
pair_coeff 2 2 -8007.8110934 0.15734738731 45.0223603063 ${R_cut}
kspace_style pppm 1.5e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes