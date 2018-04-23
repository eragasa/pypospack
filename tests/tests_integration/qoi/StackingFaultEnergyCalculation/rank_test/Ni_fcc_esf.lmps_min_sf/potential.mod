mass 1 58.6934

group Ni type 1

pair_style eam/alloy
pair_coeff * * /Users/eugeneragasa/repos/pypospack/tests/tests_integration/qoi/StackingFaultEnergyCalculation/rank_test/Ni_fcc_esf.lmps_min_sf/Ni.eam.alloy Ni

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
