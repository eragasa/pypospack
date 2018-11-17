#!/bin/bash

declare -a STRUCTURE_NAMES=(\
    "Ni_fcc_100_unit" \
    "Ni_fcc_110_unit" \
    "Ni_fcc_111_unit" \
    "Ni_bcc_100_unit" \
    "Ni_hcp_0001_unit" \
    "Ni_sc_100_unit")

for STRUCTURE_NAME in "${STRUCTURE_NAMES[@]}"
do
    IN_FILE="../${STRUCTURE_NAME}/vasp_relax_structure/CONTCAR"
    OUT_FILE="${STRUCTURE_NAME}.gga.relaxed.vasp"
    echo "$IN_FILE"
    echo "$OUT_FILE"
    python make_normalized_structures.py --in=${IN_FILE} --out=${OUT_FILE}
done
