#!/bin/bash

PYPOSPACK_SRC_DIR=/Users/eugeneragasa/repos/pypospack
API_DOC_DIR=source/pypospack/api

rm -rf $API_DOC_DIR
mkdir $API_DOC_DIR
sphinx-apidoc -f -o source/pypospack/api $PYPOSPACK_SRC_DIR
