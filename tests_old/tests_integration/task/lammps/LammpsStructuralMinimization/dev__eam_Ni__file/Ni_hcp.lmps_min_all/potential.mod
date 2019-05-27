mass 1 58.6934

group Ni type 1

pair_style eam/alloy
pair_coeff * * /home/eragasa/repos/pypospack/tests/tests_integration/task/lammps/LammpsStructuralMinimization/dev__eam_Ni__file/Ni_hcp.lmps_min_all/Ni.eam.alloy Ni

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
