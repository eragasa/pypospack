mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.507895678665762
set group O charge -1.507895678665762

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 801.8046432081591 0.32996072454302133 0.0 ${R_cut}
pair_coeff 2 2 2734.5953915915725 0.28541751651287794 28.476719144757844 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
