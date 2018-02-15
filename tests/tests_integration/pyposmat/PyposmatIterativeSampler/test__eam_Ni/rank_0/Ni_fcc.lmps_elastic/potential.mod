mass 1 58.6934

group Ni type 1

pair_style eam/alloy
pair_coeff * * /home/eugene/repos/pypospack/tests/tests_integration/pyposmat/PyposmatIterativeSampler/test__eam_Ni/rank_0/Ni_fcc.lmps_elastic/Ni.eam.alloy Ni

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
thermo	1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
