# Keppert and Wang boundary layer model

This repository contains a Fortran implementation of the Keppert and Wang
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

There is an integration test in `tests/`, run from this directory as follows:
`./tests/integration.sh`. It requires that the `namelist` has not changed
without appropriately updating the reference (`maria_WILLOUBY.*`) output in
`tests/`.

# Plan

For Maria 2017 on Puerto Rico:
- Run `footprint_model` -- done, takes about 30 CPU hours
- Run `postprocess_historical.ncl` to create a maximum wind footprint -- done,
with some moving of files

Craft a smaller example (fewer timesteps of Maria?) and create a test script
using its output -- have a version of this: `test.sh`, but it takes a couple of
minutes to complete 2 hours of Maria at dt=10s. Default timestep is dt=4s. We
could try coarsening the grid, but that's seemingly encoded in the bdy_10.txt
file (eye spatial indicies) and many of the input files -- may require running
`preprocess_historical.ncl` first.

Lint `footprint_model.f` with Fortitude

Investigate where we can tidy up
- Tidy up all the paths and organise code and data
- Improving readability (naming?)
- Factoring out functions?
- How big is the KW kernel? Should we concentrate on lifting this alone?

Possibly sack off `postprocess_historical.ncl` entirely
- Read `maria_WILLOUBY_000.d` binary files with Python and do our postprocessing there
