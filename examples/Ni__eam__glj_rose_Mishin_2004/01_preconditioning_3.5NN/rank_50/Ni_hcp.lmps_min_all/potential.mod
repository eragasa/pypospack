mass 1 58.6934

group Ni type 1

pair_style eam/alloy
pair_coeff * * /home/sullberg/repos/pypospack/examples/Ni__eam__born_exp_rose/01_preconditioning_3.5NN/rank_50/Ni_hcp.lmps_min_all/Ni.eam.alloy Ni

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
