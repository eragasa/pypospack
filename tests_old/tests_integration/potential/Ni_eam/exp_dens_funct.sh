#!/bin/bash

export PYTHONPATH=$(cd ../../;pwd)

python exp_dens_funct.py
open exp_dens_funct.png
