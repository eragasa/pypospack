mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.2688267353511766
set group O charge -2.2688267353511766

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1290.9394963687955 0.2954455082142808 0.0 ${R_cut}
pair_coeff 2 2 16513.034252320085 0.37987152256541457 62.31995581596967 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
