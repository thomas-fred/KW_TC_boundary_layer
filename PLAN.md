# Plan

[x] Version control

[x] Reproduce James' footprint plot

[x] Integration testing
- [x] 'Quick' test at 70s only a 64 core machine, but only a corner of domain
- [] Prefer a coarser grid, solved for more timesteps (therefore covering more of domain)

[x] Format / convert to modern standard
- [x] findent -- apply standard indentation, also move to free form

[] Static analysis
- [] camfort -- checks for types, implicit none completeness, dead code, etc.

[] Refactor
- [] f90split -- split program by subroutines (and other blocks?)
- [] RefactorF4Acc -- 77 into 95 with refactoring included
- [x] Multi-file build (Makefile)

[] Fix compiler warnings

[] Organise input data

[] Regenerate input data
- [] `geo_pr.nc` is on the same grid as the namelist -- how do we recreate it?
- [] Run `preprocess_historical.ncl` as part of pipeline

[] Move hardcoded parameters into namelist (e.g. reduction_factor?)

[] Improving readability (naming?)

[] Documentation?

[] Extract the kernel?

[] Sack off `postprocess_historical.ncl` entirely?
- [] Read `maria_WILLOUBY_000.d` binary files with Python and do our postprocessing there

[] Design an interface for Python

[] Locate within a workflow for coupling with other tracksets
