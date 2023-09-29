#!/bin/sh
# Usage:
#    submit.sh <Name_ID>

helpstring="Usage:
run_all.sh [Name_ID]

Name_ID: ID of the root file
"
NAME_ID=$1

# Check inputs
if [ -z ${1+x} ]
then
echo ${helpstring}
return
fi

my_string_in=$(printf "../XN/input/runID_%04d" "$NAME_ID")
my_string_out=$(printf "../XN/output/runID_%04d" "$NAME_ID")

cd ..
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
python3 modifyTree.py $my_string_in $my_string_out
