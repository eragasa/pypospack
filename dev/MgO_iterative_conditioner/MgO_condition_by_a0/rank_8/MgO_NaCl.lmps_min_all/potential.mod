mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.6287378290743892
set group O charge -1.6287378290743892

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1124.7578630014054 0.3100029336244492 0.0 ${R_cut}
pair_coeff 2 2 19494.628108092067 0.39944549942886554 76.09156141505969 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
