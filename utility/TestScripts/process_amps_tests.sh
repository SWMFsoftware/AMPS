#!/bin/csh
cd $HOME/Sites

# Remove code from yesterday as well as various logs
rm -rf test.diff

# Create directory for new test results
setenv NEWTESTDIR AMPS_TEST_RESULTS/`date -v-1d +%Y/%m/%d`
mkdir -p ${NEWTESTDIR}

# Copy over test results but preserve the subdirectories in Current
chmod a+r Current/*/*
cp -r Current/* ${NEWTESTDIR}/
rm -f Current/*/*

# Create new index.html and list changes in test.diff
./process_amps_tests.pl
