#!/bin/bash

# we have ifort on ARC
# ifort -O3 -cpp -qopenmp -ftz -zero -o footprint_model footprint_model.f

# but not on pwlldu, so use GNU fortran
# stds '2008' and '2018' both fail with errors, so using 'gnu'
gfortran -O3 -cpp -std='gnu' -fopenmp -o footprint_model footprint_model.f
