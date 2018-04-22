mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.3722212935477476
set group O charge -2.3722212935477476

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1015.3903115974379 0.30888210713117226 0.0 ${R_cut}
pair_coeff 2 2 4225.977820318988 0.21074584996471502 76.13878777141537 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
