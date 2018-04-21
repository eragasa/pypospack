mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 2.0754071601519097
set group O charge -2.0754071601519097

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 916.7503018596776 0.31019793735679424 0.0 ${R_cut}
pair_coeff 2 2 1585.568223002891 0.3273790901546082 60.424639218392905 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
