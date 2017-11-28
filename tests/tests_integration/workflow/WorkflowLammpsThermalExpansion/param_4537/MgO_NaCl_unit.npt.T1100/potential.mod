mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.68595133453
set group O charge -1.68595133453

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1257.74769115 0.28905643469 0.0 ${R_cut}
pair_coeff 2 2 3064.42109336 0.236064465148 18.7313535189 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
