mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.64252902328
set group O charge -1.64252902328

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1425.65456865 0.277130619093 0.0 ${R_cut}
pair_coeff 2 2 8710.58674493 0.211543782712 -4.62088168449 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
