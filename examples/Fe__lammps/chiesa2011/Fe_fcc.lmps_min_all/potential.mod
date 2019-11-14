mass 1 55.845

group Fe type 1

pair_style eam/alloy
pair_coeff * * chiesa2011.alloy Fe

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
