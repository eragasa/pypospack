mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.63539617609
set group O charge -1.63539617609

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1445.86118726 0.28350177361 0.0 ${R_cut}
pair_coeff 2 2 14001.7361193 0.0970682303871 13.3227844282 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
