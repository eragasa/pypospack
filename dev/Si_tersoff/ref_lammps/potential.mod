mass 1 28.0855

group silicon type 1

variable R_cut equal 10.0

pair_style      tersoff/mod
pair_coeff 	* * Si.tersoff.mod Si 
#pair_coeff 	* * Si.tersoff.mod Si Si Si Si Si Si Si Si

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
