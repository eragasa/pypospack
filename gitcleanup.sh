#!/bin/bash
find . | grep -E "(__pycache__|\.pyc|\.pyo$)" | xargs rm -rf
#find . | grep -E "(\.lmps_elastic|\.lmps_min_all)" | xargs rm -rf
