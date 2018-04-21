mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.2951801394729396
set group O charge -2.2951801394729396

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 883.4667805417473 0.31166982185035813 0.0 ${R_cut}
pair_coeff 2 2 19911.75479118061 0.2926020667790446 42.879084588207746 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
