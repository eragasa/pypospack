mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.6421645765040513
set group O charge -1.6421645765040513

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 976.3962522474803 0.3291762116597489 0.0 ${R_cut}
pair_coeff 2 2 3230.583724315258 0.1848642364786558 56.128495319422115 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
