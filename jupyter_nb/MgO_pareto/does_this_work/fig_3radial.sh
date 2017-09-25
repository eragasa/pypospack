#!/bin/bash
export PYTHONPATH=$(cd ../../pyflamestk/;pwd)
python fig_3radial.py
open fig_3radar.png
