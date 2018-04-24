mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.0306397248688364
set group O charge -2.0306397248688364

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 984.9627421449446 0.2926162890476923 0.0 ${R_cut}
pair_coeff 2 2 9364.731711034323 0.11177352626295374 67.67542120127403 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
