#!/bin/bash

# we have ifort on ARC
# ifort -O3 -cpp -qopenmp -ftz -zero -o footprint_model footprint_model.f

# but not on pwlldu, so use GNU fortran
rm footprint_model
gfortran -O3 -cpp -std='gnu' -fopenmp -o footprint_model src/main.f
