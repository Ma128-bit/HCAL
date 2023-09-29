#!/bin/sh
# Usage:
#    hadd_run.sh <Name>

helpstring="Usage:
hadd_run.sh [Name]\n
Name: Name of the root file with no .root
"
NAME=$1

# Check inputs
if [ -z ${1+x} ]
then
echo ${helpstring}
return
fi
cmsenv
cd $NAME
hadd "${NAME}_new.root" ./output/*
cd ..
