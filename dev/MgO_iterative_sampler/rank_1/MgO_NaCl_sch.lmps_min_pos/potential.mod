mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.6102951165102182
set group O charge -1.6102951165102182

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 872.8291467350588 0.29889716378497794 0.0 ${R_cut}
pair_coeff 2 2 3727.0736816732065 0.3476815772257621 63.9105328410665 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
