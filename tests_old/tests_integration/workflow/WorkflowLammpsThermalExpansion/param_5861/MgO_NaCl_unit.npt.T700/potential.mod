mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.63577302612
set group O charge -1.63577302612

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1558.72850733 0.278190610155 0.0 ${R_cut}
pair_coeff 2 2 26609.5495557 0.148782776904 9.00084075824 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
