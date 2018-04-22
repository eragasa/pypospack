mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.6733610086346575
set group O charge -1.6733610086346575

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1185.910957161503 0.30936047509394593 0.0 ${R_cut}
pair_coeff 2 2 19776.92883027598 0.15101676301755923 25.910104919364944 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
