#!/bin/bash
export PYTHONPATH=$(cd ../../pyflamestk/;pwd)
python fig_combo.py
open fig_combo.png
