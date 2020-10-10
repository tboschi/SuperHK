#! /bin/bash

usage="Usage: $0 info_file

Check if jobs were completed and repair if needed
Works with HTCondor or Slurm.

  parameters
    info_file    path to file with .info extension in output folder"

if [ "$#" -ne 1 ] ; then
	echo This script requires one argument. Check usage with
	echo $0 -h
	exit 1
elif [ "$1" == "-h" ] ; then
	echo "$usage" >&2
elif [ ! -s "$1" ] ; then
	echo Argument \"$1\" is not valid. Check usage with
	echo $0 -h
	exit 1
fi

SCHED=""
if condor_q &> /dev/null ; then
	SCHED="HTCONDOR"
	sub=condor_submit
elif squeue &> /dev/null ; then
	SCHED="SLURM"
	sub=sbatch
else
	echo There is neither HTCondor nor Slurm on this machine. I am sorry, I cannot help you
	exit 1
fi


Sens=$PWD/cross-fitter.sh

name=${1%.*}

while read -r point ; do
	# if finished is not shown, means something bad happened
	echo Checking point $point

	outdir=$name'_'$point
	outdir=$(realpath $outdir)

	allog=$(find $outdir -name "*.log")
	allog=(${allog})
	all=${#allog[@]}

	if [ "$all" -eq 0 ]; then
		echo No log files for point $point, job did not even start. Restarting
		if [ -s "$outdir/R"*".$point.sub" ] ; then
			$sub $outdir/R*.$point.sub
		fi
	else
		for log in "${allog[@]}" ; do
			if ! tail $log | grep -q Finished ; then
				#log is #/long/path/to/folder/CPV_scan_point/fitter.proc.log
				proc=${log##*/}
				#proc is fitter.proc.log
				proc=${proc#*.}
				#proc is proc.log
				proc=${proc%.*}
				#proc is proc
				echo process $point : $proc did not complete

				scriptname=$outdir/Recover.$proc.sub
				if [ "$SCHED" == "HTCONDOR" ] ; then
					cat > $scriptname << EOF
#! /bin/bash
# script submission for SLURM
# sumbit with --
#	$sub $scriptname

executable		= $Sens
arguments		= $proc $all $outdir/this_sensitivity.card
getenv			= True
#requirements		= HasAVXInstructions
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $log
error			= $log

queue
EOF
				elif [ "$SCHED" == "SLURM" ] ; then
					cat > $scriptname << EOF
#! /bin/bash
# script submission for Slurm
# sumbit with --
#	$sub $scriptname

#SBATCH --job-name=fixing
#SBATCH -o $log
#SBATCH -p nms_research,shared
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1

srun $Sens $proc $all $outdir/this_sensitivity.card

EOF
				fi
				$sub $scriptname
			fi
		done
	fi
done < $1
