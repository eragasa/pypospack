mass 1 58.6934

group Ni type 1

pair_style eam/alloy
pair_coeff * * /ufrc/phillpot/eragasa/pareto_optimization_3.5NN/rank_21/Ni_fcc_110_s.lmps_min_pos/Ni.eam.alloy Ni

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes