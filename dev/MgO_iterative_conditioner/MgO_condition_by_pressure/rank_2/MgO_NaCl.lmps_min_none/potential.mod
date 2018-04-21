mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.9960342126287767
set group O charge -1.9960342126287767

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1073.1341673153372 0.3280360093167101 0.0 ${R_cut}
pair_coeff 2 2 20867.59845877755 0.2998654970783864 41.122990685443455 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
