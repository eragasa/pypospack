#!/bin/bash
if [ -e "out.dat" ]
then
  rm out.dat
fi

THIS_HOSTNAME=$(hostname -f)
if [ "$THIS_HOSTNAME" = minerva ]; then
   # settings for Eugene's laptop
   export LAMMPS_BIN=/usr/local/bin/lammps
fi

# run lammps simulation
${LAMMPS_BIN} -i in.min_pos > out.dat
