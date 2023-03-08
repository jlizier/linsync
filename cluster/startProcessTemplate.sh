#!/bin/bash
#
# Linear Sync Toolkit (linsync)
# Copyright (C) 2023 Joseph T. Lizier
# Distributed under GNU General Public License v3
# Modified from Simple cluster jobs framework distributed under GPL by J.T. Lizier
#
# -----------------------------------------------------------------
#
# Main script to submit to the cluster.
# Sets up some resources and environment variables and then starts the
#  Matlab script.
#
# For some clusters this sets the current working directory to here:
#$ -cwd

# 0. Set the PBS properties: (These are set for my Usyd cluster)
#PBS -P MultiVInt
#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=04:00:00
#PBS -M joseph.lizier@sydney.edu.au

# 1. Set your current working directory for the process.
# For some clusters this is done with "$ -cwd" above.
# For job submission in other cluster types, use the following line:
cd $PBS_O_WORKDIR

# 2. Load any modules or set any paths you require.
module add matlab/R2019a

# 3. Template for starting your process:
# You can use [@PROPSFILE], [@LOGFILE] and [@ERRFILE] to refer to the specific
#  input properties file and log files created/assigned for this process.
echo "runComputeSyncResults('[@PROPSFILE]');" | matlab -nodesktop > [@LOGFILE] 2> [@ERRFILE]

