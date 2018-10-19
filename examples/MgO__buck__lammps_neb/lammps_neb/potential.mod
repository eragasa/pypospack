mass 1 24.305
mass 2 15.999

group Mg type 1
group O type 2

set group Mg charge 1.66039996869
set group O charge -1.66039996869

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1528.2261267 0.282600919639 0.0 ${R_cut}
pair_coeff 2 2 6677.06285722 0.138925130231 8.13419383626 ${R_cut}

kspace_style pppm 1.5e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes