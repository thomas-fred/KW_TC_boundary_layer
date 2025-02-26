# Keppert and Wang boundary layer model

This repository contains a Fortran implementation of the Keppert and Wang
(2001) boundary layer model for coupling gradient-level tropical cyclone winds
to the earth's surface.

The boundary layer model was designed and coded by Yuqing Wang.
The code was modified by James Done and Ming Ge in 2019-2021.

# Aims

Run this model for intense synthetic tropical cyclones for Jamaica or Mauritius

# Usage

The code is currently split into three files.

## preprocess_historical.ncl

This program selects, processes, and plots all TC tracks within 200km of land
Input files:  rmax_ib.nc, vmax_ib.nc, IBTrACS.ALL.v04r00.nc, geo_pr.nc
Output files: bdy_10min.txt, lat_2d.dat, topo.dat, landuse.dat

In the bundle from James, these outputs are already given for Puerto Rico, so
we don't need to run preprocess_historical.ncl

## footprint_model

I'm using gfortran, James uses ifort. Compile with `./compile.sh`

Set parameters in `namelist`. Note that many parameters are still hardcoded
inside `fooprint_model.f`, e.g. filepaths.

Run with `./footprint_model`

This will run for the timesteps requested in `namelist` and will produce lots
of binary files like `maria_WILLOUBY_000.d`

## postprocess_historical.ncl

TODO

# Plan

For Maria 2017 on Puerto Rico:
- Run `footprint_model` -- done, takes about 30 CPU hours
- Run `postprocess_historical.ncl` to create a maximum wind footprint

Craft a smaller example (fewer timesteps of Maria?) and create a test script
using its output

Lint `footprint_model.f` with Fortitude

Investigate where we can tidy up
- Tidy up all the paths and organise code and data
- Improving readability (naming?)
- Factoring out functions?
- How big is the KW kernel? Should we concentrate on lifting this alone?

Possibly sack off `postprocess_historical.ncl` entirely
- Read `maria_WILLOUBY_000.d` binary files with Python and do our postprocessing there
