#!/bin/ksh
#
# Linear Sync Toolkit (linsync)
# Copyright (C) 2023 Joseph T. Lizier
# Distributed under GNU General Public License v3
# Originally distributed as Simple cluster jobs framework under GPL by J.T. Lizier
#
# -----------------------------------------------------------------
#
# delete all of our submitted/running jobs, or those with sched ids above $1
#

schedIds=$(qstat | grep $(whoami) | sed 's/^\([0-9]*\)\..*$/\1/')

for schedId in $schedIds; do
	if (($# > 0)); then
		# we're only qdeling for sched id's above $1
		if  (($schedId >= $1)); then
			echo "qdel-ing $schedId"
			qdel $schedId
		else
			echo "not qdel-ing $schedId"
		fi
	else
		echo "qdel-ing $schedId"
		qdel $schedId
	fi
done

