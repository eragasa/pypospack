mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.3110959979317283
set group O charge -2.3110959979317283

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1017.897159500457 0.30562958336551893 0.0 ${R_cut}
pair_coeff 2 2 8626.82442584353 0.3267719247097669 27.199594411545267 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
