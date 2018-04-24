mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.5539689568370378
set group O charge -1.5539689568370378

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1251.111218423101 0.2904898324677775 0.0 ${R_cut}
pair_coeff 2 2 17791.26594453386 0.3445982013966401 54.369572613819685 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
