#!/bin/bash
#
# Linear Sync Toolkit (linsync)
# Copyright (C) 2023 Joseph T. Lizier
# Distributed under GNU General Public License v3
# Modified from Simple cluster jobs framework distributed under GPL by J.T. Lizier
#
# -----------------------------------------------------------------
#
# Submits one job of startProcessTemplate.sh to the cluster, 
#  which calls runComputeSyncResults.m with parameters.m
#  fed as the parameters file after substituting command line arguments
#  $1, $2 and $3 for [@P1], [@P2] and [@P3] in the text respectively.
#
# Args:
# $1 substitute for [@P1] in parameters.m (optional)
# $2 substitute for [@P2] in parameters.m (optional)
# $3 substitute for [@P3] in parameters.m (optional)


datesec=$(date +%s)
datenano=$(date +%N)
# Make the run id the current time in millisec
runId=j${datesec}${datenano:0:3}
# Set up log files and temp run files using this job id
logLocalname="$runId".log
logfile=log/$logLocalname
logfileEscaped=log\\/$logLocalname
errLocalname="$runId".err
errfile=err/$errLocalname
errfileEscaped=err\\/$errLocalname
propsLocalname=properties"$runId".m
propsfile=props/$propsLocalname
propsfileEscaped=props\\/$propsLocalname
processLocalname="$runId".sh
processfile=processshellscripts/$processLocalname

# Create properties file
echo "substituting $1, $2 and $3"
sed "s|\[@P1\]|$1|g" parameters.m > genInput.properties
sed "s|\[@P2\]|$2|g" genInput.properties > genInput2.properties
sed "s|\[@P3\]|$3|g" genInput2.properties > $propsfile
rm genInput2.properties
rm genInput.properties

sed "s/\[@LOGFILE\]/$logfileEscaped/g" startProcessTemplate.sh > startProcess.sh
sed "s/\[@ERRFILE\]/$errfileEscaped/g" startProcess.sh > startProcess2.sh
sed "s/\[@PROPSFILE\]/$propsfileEscaped/g" startProcess2.sh > $processfile
rm startProcess2.sh
rm startProcess.sh

chmod u+x $processfile

# Submit the job:
echo
echo "submitting $processfile ..."

######
# YOU WILL NEED TO SELECT THE APPROPRIATE LINE DEPENDING ON YOUR SYSTEM HERE
#  AND EDIT THE RUNTIME AND MEMORY REQUIREMENTS IF REQUIRED
###### 
# On some cluster types specify resources in qsub:
#qsub -lvmem=500MB,walltime=23:59:00 $processfile
# On other cluster types: (default on my system was 24 hours runtime - change with -l long (get 7 days, but lower priority)
# qsub -cwd -l mf=3.5G $processfile
# Or (USyd cluster) no resources, they're specified in the script separately:
qsub $processfile

echo
echo "To tail the output ...:"
echo
echo "tail -f $logfile"

# sleep 2
# tail -f $logfile

echo off


