mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.533842542173621
set group O charge -1.533842542173621

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 952.0284891774219 0.3067179239266943 0.0 ${R_cut}
pair_coeff 2 2 24516.907297092912 0.15975517110357593 41.68420575310288 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
