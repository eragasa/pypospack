mass 1 26.982

group Al type 1

pair_style eam/alloy
pair_coeff * * /home/sullberg/repos/pypospack/examples/Al__eam__born_exp_fs/pareto_optimization_3.5NN/rank_2/Al_fcc_100_s.lmps_min_pos/Al.eam.alloy Al

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
