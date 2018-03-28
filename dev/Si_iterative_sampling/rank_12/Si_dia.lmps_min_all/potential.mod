mass 1 28.086

group Si type 1

pair_style sw
pair_coeff * * Si.parameters Si


neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
