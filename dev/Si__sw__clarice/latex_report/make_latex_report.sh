#!/bin/bash

python make_latex_report.py
latexmk -pdf -pv latex_report_test

