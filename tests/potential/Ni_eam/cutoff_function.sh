#!/bin/bash

export PYTHONPATH=$(cd ../../;pwd)
python cutoff_function.py
open cutoff_function.png
