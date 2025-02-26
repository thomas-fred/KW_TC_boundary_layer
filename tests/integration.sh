#!/bin/bash

# Test footprint_model and postprocess_historical.ncl against reference output
# To be run from directory containing `tests/` as follows:
# ./tests/integration.sh

set -e

echo "Compiling program"
./compile.sh

echo "Using following problem configuration"
cat namelist

STORM_YEAR="MARIA_2017"
rm -r ${STORM_YEAR}
mkdir ${STORM_YEAR}

echo "Running boundary layer simulation"
./footprint_model
# the postprocessing script will try and work on all the rows in the bdy_10min.txt file
# our namelist is specifying a shorter simulation, so truncate the bdy_10min.txt first
# number of lines = simulation hours * 6 + 1 (header) + 1 (for luck?)
KFI_KEY_VALUE=$(cat namelist | grep kfi)
DURATION_HOURS=${KFI_KEY_VALUE#*=}
N_ROWS_TRACK_FILE=$(echo "${DURATION_HOURS} * 6 + 2" | bc)
head -${N_ROWS_TRACK_FILE} bdy_10min.txt > ${STORM_YEAR}/bdy_10min.txt
mv *_???.d ${STORM_YEAR}

echo "Postprocessing outputs; creating footprint plot"
micromamba run --name kw-pbl ncl postprocess_historical.ncl

echo "Checking output against reference"
STORM_PROFILE="maria_WILLOUBY"
# with set -e, any non-zero exit code should terminate early
diff ${STORM_YEAR}/${STORM_PROFILE}.txt tests/${STORM_PROFILE}.txt
./tests/visual_compare.sh ${STORM_YEAR}/${STORM_PROFILE}.pdf tests/${STORM_PROFILE}.pdf

echo "Passed!"
