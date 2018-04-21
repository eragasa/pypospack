mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.3002490590750533
set group O charge -2.3002490590750533

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1247.2875025704816 0.3251901441223064 0.0 ${R_cut}
pair_coeff 2 2 11985.964246506135 0.26471397353488824 68.11341769118575 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
