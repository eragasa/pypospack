mass 1 58.6934

group Ni type 1

pair_style eam/alloy
pair_coeff * * /home/sullberg/repos/pypospack/dev/cluster_sampling/dev__cluster_sampling/rank_5/Ni_fcc_vac.lmps_min_pos/Ni.eam.alloy Ni

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
