mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.673390077367403
set group O charge -1.673390077367403

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 991.1403258979893 0.3226693822023098 0.0 ${R_cut}
pair_coeff 2 2 6391.819978791292 0.3643994700619103 75.3455597303464 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
