mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.9179406879711596
set group O charge -1.9179406879711596

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 810.4136342137181 0.29896087357084233 0.0 ${R_cut}
pair_coeff 2 2 4888.921747208894 0.3126294278687661 51.00665792627922 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
