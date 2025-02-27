# Kepert and Wang boundary layer model

This repository contains a Fortran implementation of the Kepert and Wang
(2001) boundary layer model for coupling gradient-level tropical cyclone winds
to the earth's surface.

The boundary layer model was designed and coded by Yuqing Wang.
The code was modified by James Done and Ming Ge in 2019-2021.

# Aims

Improve the usability of this model
Run this model for intense synthetic tropical cyclones for Jamaica or Mauritius

# Usage

The code is currently split into three files. The first and last are for pre
and post processing and are written in NCAR Command Language (NCL), a scripting
language developed at NCAR and now not maintained. The other contains the
boundary layer model, written in Fortran of an unknown standard.

To install NCL and other dependencies:
```
micromamba create -f environment.yml -y
```

To make these available:
```
micromamba activate kw-pbl
```

You will also need a fortran compiler, e.g. gfortran or ifort.

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

# Testing

There are simple integration tests in `tests/integration`, run from this directory as follows.

For 2 hours of Maria 2017 over Puerto Rico, with a coarse timestep:
`./tests/integration/test.sh maria_PRI_short`

The same problem, but 24 hours and with a 4 second timestep:
`./tests/integration/test.sh maria_PRI_long`
