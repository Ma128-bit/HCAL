# MPGD-HCAL signal fit
This repository contains the files necessary to perform the fit on the signals on each strip/pad with Condor.

## Setting the environment

```
cmsrel CMSSW_13_2_4
cd CMSSW_13_2_4/src
cmsenv
git clone https://github.com/Ma128-bit/HCAL.git .
source setup.sh
```
<p>&nbsp;</p>

In the `Root_File` directory you need to copy the files you want to analyze

## Run the analysis 
To analyze a file:
```
source prepare_run.sh [file_name] [events_per_file]
```
* [file_name] = name of the file without .root
* [events_per_file] = number ov events per subfile

This crate a directory named [file_name] with the following items:
* `input` a directory containing the subfiles into which the initial file was divided (each consisting of [events_per_file] events)
* `output` and `log` empty directories
* `launch_analysis.sh` a file that launch on condor the analysis of one subfile
  * Example: `source launch_analysis.sh 5` submit the analysis of the subfile number 5
* `submit.condor` used to submit the analysis of all the subfiles
