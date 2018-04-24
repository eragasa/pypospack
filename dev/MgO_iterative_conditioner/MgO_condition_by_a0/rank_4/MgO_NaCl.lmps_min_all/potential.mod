mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.8013157705069975
set group O charge -1.8013157705069975

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 800.2098472747266 0.3033662810872803 0.0 ${R_cut}
pair_coeff 2 2 14632.037663109702 0.23400317648039198 69.66619928259445 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
