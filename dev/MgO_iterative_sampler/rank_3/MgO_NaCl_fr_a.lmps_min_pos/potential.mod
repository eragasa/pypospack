mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.713768000698136
set group O charge -1.713768000698136

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1027.6516591691873 0.2934831246256318 0.0 ${R_cut}
pair_coeff 2 2 21120.842859678647 0.18037284948043203 60.17316407089181 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
