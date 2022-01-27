#!/bin/bash
#Wrapper developed by Dr. Sreenu Vattipally at University of Glasgow, adapted by Simon Daldry at University of Glasgow

cores="7";

usage() { echo "Usage: $0 -i inputDirectory -c cores" 1>&2; exit 1; }

if [ "$#" -lt 1 ]; then
usage;
fi
while getopts ":i:c:" opt; do
    case "${opt}" in
        i) runDir=${OPTARG} ;;
        c) cores=${OPTARG} ;;
        *) usage ;;
    esac
done
shift "$((OPTIND-1))"

ulimit -n 4000

while true
do
	runStatus=$(grep CompletionStatus $runDir/RunCompletionStatus.xml |grep -c CompletedAsPlanned )
	if [ $runStatus -eq 1 ]; then
		snakemake -j $cores
		exit;
	else
		echo "Waiting for sequencing to complete...";
		sleep 180;			# Wait for three minutes
	fi
done
