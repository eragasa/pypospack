#!/bin/bash
rm out.dat
source /opt/intel/composer_xe_2011_sp1.8.273/bin/compilervars.sh 'intel64'
/home/eugene/bin/lmp_intx_comb -i in.elastic  > out.dat
