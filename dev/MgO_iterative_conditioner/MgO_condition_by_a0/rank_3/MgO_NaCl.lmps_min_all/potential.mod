mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.175312172413404
set group O charge -2.175312172413404

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 882.0519485284365 0.29258186603958153 0.0 ${R_cut}
pair_coeff 2 2 11003.060129170797 0.23570520996424668 57.21963420129103 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
