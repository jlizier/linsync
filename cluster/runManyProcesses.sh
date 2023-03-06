#!/bin/ksh
#
# Linear Sync Toolkit (linsync)
# Copyright (C) 2023 Joseph T. Lizier
# Distributed under GNU General Public License v3
# Modified from Simple cluster jobs framework distributed under GPL by J.T. Lizier
#
# -----------------------------------------------------------------
#
# Example script to submit several processes with one parameter [@P1] updated
#  in the parameters.m file for each run
#
# Args:
# $1 start parameter [@P1]
# $2 increment
# $3 end parameter

start=$1
increment=$2
end=$3

current=$start
while ((current <= end)); do
	echo "Running for $current"
	./runProcess.sh $current
	# Ensure we sleep for >1ms to increment job number properly:
	sleep 0.002
	current=$((current + increment))
done

# Alternate coding for this for other shell:
#while [[ $(echo "$current <= $end" | bc) -ne 0 ]]; do
#	echo "$current $increment"
#	./runRbns.sh $current
#	# Ensure we sleep for 1ms to increment job number properly:
#	sleep 0.002
#	current=$(echo "$current+$increment" | bc)
#done

sleep 0.002
