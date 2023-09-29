#!/bin/sh
# Usage:
#    prepare_run.sh <Name_file> <N_files>

helpstring="Usage:
prepare_run.sh [Name_file] [N_files]
\n-Name_file: Name of the initaial root file without .root
\n-N_files: number of files into which to split the initial file
"
NAME=$1
NUMBER=$2

# Check inputs
if [ -z ${2+x} ]
then
echo -e ${helpstring}
return
fi

cmsenv
./split ${NAME}.root ${NUMBER}

directory="./$NAME/input/"
#num_files=`ls -l $directory | grep -v ^d | wc -l`
num_files=$(find "$directory" -maxdepth 1 -type f | wc -l)
#echo "Number of files: $num_files"

cp ./templates/submit.condor ./$NAME
echo "queue $num_files" >> "./$NAME/submit.condor"
chmod a+x ./$NAME/submit.condor

cp ./templates/launch_analysis.sh ./$NAME
sed -i "s/XN/$NAME/g" ./$NAME/launch_analysis.sh
chmod a+x ./$NAME/launch_analysis.sh

