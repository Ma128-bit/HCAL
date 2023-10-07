# MPGD-HCAL signal fit
This repository contains the files necessary to perform the fit on the signals on each strip/pad with Condor.

## Setting the environment
**NOTE**: The use of CMSSW is **not** mandatory, mainly ROOT and python3 are required, however, if CMSSW is not used: 
* `cmsenv` must be removed from hadd_run.sh, prepare_run.sh, setup.sh and templates/launch_analysis.sh
* `source /cvmfs/cms.cern.ch/cmsset_default.sh` must be removed from templates/launch_analysis.sh
* Probably other changes are needed based on where the code is running.

These instructions assume the use of CMSSW:
```
cmsrel CMSSW_13_2_4
cd CMSSW_13_2_4/src
cmsenv
git clone https://github.com/Ma128-bit/HCAL.git .
source setup.sh
```
In the `Root_File` directory you need to copy the files you want to analyze.
<p>&nbsp;</p>

## Run the analysis 
To analyze a file:
```
source prepare_run.sh [file_name] [events_per_file]
```
* [file_name] = name of the file without .root
* [events_per_file] = number of events per subfile

This crate a directory named [file_name] with the following items:
* `input` a directory containing the subfiles into which the initial file was divided (each consisting of [events_per_file] events)
* `output` and `log` empty directories
* `launch_analysis.sh` a file that launch on condor the analysis of one subfile
  * Example: `source launch_analysis.sh 5` submit the analysis of the subfile number 5
* `submit.condor` used to submit the analysis of all the subfiles

<p>&nbsp;</p>

## After completing the analysis
In the `[file_name]/output` directory you will find the modified subfiles and in `[file_name]/log` the logs associated with the submission to condor.
To add all the output subfiles:
```
source hadd_run.sh.sh [file_name]
```
This will produce the final file in the [file_name] directory.

<p>&nbsp;</p>
