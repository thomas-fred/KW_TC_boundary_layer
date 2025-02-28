#!/bin/bash

# Test boundary_layer and postprocess_historical.ncl against reference output
# To be run from directory containing `tests/`, e.g.
# ./tests/integration/test.sh maria_PRI_short

TEST_START_TIME=$EPOCHREALTIME

TEST_DIR="tests/integration/$1"
if [ ! -d ${TEST_DIR} ]; then
    echo "${TEST_DIR} does not exist, quitting"
    exit 1
fi

STORM_YEAR="MARIA_2017"
rm -r ${STORM_YEAR}
mkdir ${STORM_YEAR}

echo "Compiling program"
cd src
make clean
make
make install
cd ..

echo "Using following problem configuration for ${STORM_YEAR}"
TEMP_NAMELIST_FILENAME=$(mktemp)
cat ./namelist > ${TEMP_NAMELIST_FILENAME}
rm ./namelist
cp "${TEST_DIR}/namelist" "./namelist"
cat ./namelist

echo "Running boundary layer simulation"
./boundary_layer

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

# keep track of whether we're in good shape or not
TEST_STATE=0

diff ${STORM_YEAR}/${STORM_PROFILE}.txt ${TEST_DIR}/${STORM_PROFILE}.txt &>/dev/null
if [ $? != 0 ]; then
    echo "Failure: mismatch in computed maximum speeds"
    TEST_STATE=1
fi

./tests/visual_compare.sh ${STORM_YEAR}/${STORM_PROFILE}.pdf ${TEST_DIR}/${STORM_PROFILE}.pdf
if [ $? != 0 ]; then
    echo "Failure: mismatch in generated plot"
    TEST_STATE=1
fi

# replace original namelist
rm ./namelist
cp ${TEMP_NAMELIST_FILENAME} ./namelist
rm ${TEMP_NAMELIST_FILENAME}

ELAPSED_TIME=$(bc -l <<< "$EPOCHREALTIME - $TEST_START_TIME")
echo "Test took ${ELAPSED_TIME} seconds"

if [ $TEST_STATE != 0 ]; then
    echo "Test failed :("
    exit $TEST_STATE
else
    echo "Test passed :)"
    exit $TEST_STATE
fi
