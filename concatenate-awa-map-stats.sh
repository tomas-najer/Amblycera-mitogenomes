#!/bin/bash

sequence=$1 # SRR number
subsampling=$2 # 2mil or 10mil

function part4 {
	for k in 0 5 10 15 20 25 30 35 40; do
		echo "$subsampling-$sequence-${k}-mer.log" >> $subsampling-$sequence-final.map.log
		tail -n 7 $subsampling-$sequence-${k}-mer.log >> $subsampling-$sequence-final.map.log
	done
}

part4
exit

# save as concatenate-map-stats.sh, run chmod 755 *sh and run the script with sequence number and subsampling
