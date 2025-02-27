#!/bin/bash

# Test footprint_model and postprocess_historical.ncl against reference output
# To be run from directory containing `tests/`, e.g.
# ./tests/integration/test.sh maria_PRI_short

TEST_DIR="tests/integration/$1"
if [ ! -d ${TEST_DIR} ]; then
    echo "${TEST_DIR} does not exist, quitting"
    exit 1
fi

echo "Compiling program"
TEMP_NAMELIST_FILENAME=$(mktemp)
cat ./namelist > ${TEMP_NAMELIST_FILENAME}
rm ./namelist
cp "${TEST_DIR}/namelist" "./namelist"
./compile.sh

echo "Using following problem configuration"
cat ./namelist

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
micromamba run --name kw-pbl ncl ./scripts/postprocess_historical.ncl

echo "Checking output against reference"
STORM_PROFILE="maria_WILLOUBY"

diff ${STORM_YEAR}/${STORM_PROFILE}.txt ${TEST_DIR}/${STORM_PROFILE}.txt &>/dev/null
test $? -eq 0 || echo "Failure: mismatch in computed maximum speeds"

./tests/visual_compare.sh ${STORM_YEAR}/${STORM_PROFILE}.pdf ${TEST_DIR}/${STORM_PROFILE}.pdf
test $? -eq 0 || echo "Failure: mismatch in generated plot"

# replace original namelist
rm ./namelist
cp ${TEMP_NAMELIST_FILENAME} ./namelist
rm ${TEMP_NAMELIST_FILENAME}
