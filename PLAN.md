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

Lint
- camfort -- checks for types and dead code

Format
- ffixed2free -- should move us much closer to a more modern standard like f95
- fortify -- apply standard indentation

Refactor
- f90split -- split program by subroutines (and other blocks?)
- RefactorF4Acc -- 77 into 95 with refactoring included

Investigate where we can tidy up
- Tidy up all the paths and organise code and data
- Improving readability (naming?)
- Factoring out functions?
- How big is the KW kernel? Should we concentrate on lifting this alone?

Possibly sack off `postprocess_historical.ncl` entirely
- Read `maria_WILLOUBY_000.d` binary files with Python and do our postprocessing there
