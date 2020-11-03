! /bin/bash

usage="Usage: $0 info_file [-f]

Check if jobs were completed and repair if needed
Works with HTCondor or Slurm.

  parameters
    info_file    path to file with .info extension in output folder
    -f		 repair files"

rep=false;
while getopts 'fh' flag; do
	case "${flag}" in
		f) rep=true;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done

shift $((OPTIND-1))

if [ "$#" -ne 1 ] ; then
	echo This script requires one argument. Check usage with
	echo $0 -h
	exit 1
elif [ ! -s "$1" ] ; then
	echo Argument \"$1\" is not valid. Check usage with
	echo $0 -h
	exit 1
fi


Sens=$PWD/cross-fitter.sh
nameExec=${Sens##*/}
nameExec=${nameExec%.*}
trisens=trisens.sh


if pgrep $trisens ; then
	echo Trisens is still running. Wait for it to finish.
	exit 1
fi


SCHED=""
queue=".queue_status"
if condor_q &> /dev/null ; then
	SCHED="HTCONDOR"
	sub=condor_submit
	condor_q -autoformat cmd args > $queue
elif squeue &> /dev/null ; then
	SCHED="SLURM"
	sub=sbatch
	squeue -h -r -u $USER -o \"%o %K\" > $queue
else
	echo There is neither HTCondor nor Slurm on this machine. I am sorry, I cannot help you
	exit 1
fi

name=${1%.*}
root=${nam%/*}

repeat=()
while read -r point ; do
	# if finished is not shown, means something bad happened
	echo Checking point $point

	outdir=$name'_'$point
	outdir=$(realpath $outdir)
	script=$outdir/R$nameExec.$point.sub

	# directory does not exists or no script
	if ! [ -d $outdir ] || ! [ -s $script ] ;
		echo Detected: directory $outdir or $script do not exist
		if [ "$rep" == "true" ] ; then
			echo Point $point will be resubmitted
			repeat=(${repeat[@]} $point)
		else
			echo Point $point is missing. You can resubmit with -f to repair
		fi
		continue
	fi

	# find all root files
	out=$(find $outdir -name "*.root")
	out=(${out})

	if [ ${#out[@]} -eq 0 ] ; then	# no root files
		# check if job is still running
		echo Detected: there are no output files in $outdir
		if grep -q $point $queue; then 
			echo $nameExec for point $point is still running. Just wait.
		elif [ "$rep" == "true" ] ; then
			echo Point $point will be resubmitted
			repeat=(${repeat[@]} $point)
		else
			echo $nameExec for point $point failed and never ran. You can resubmit with -f to repair
		fi
		continue
	fi 


	# find how many jobs there should be
	if grep -q queue $script ; then
		job=$(grep queue $script | cut -f2 -d' ')
	elif grep -q array $script ; then
		job=$(grep array $script | cut -f4 -d'-')
		job=$((job + 1))	# arrays start from 0
	fi

	repair=()

	# now go through each file individually
	for ff in "${out[@]}" ; do
		# this is the file number
		num=${num%.root}
		num=${num##*.}

		log=$outdir/L$nameExec.$num.log

		# remove wrong stuff
		if [ $num -ge $job ] ; then
			rm -f $ff $log
			continue
		fi

		# if root file is older than script
		if [ $ff -ot $script ] ; then	# maybe failed?
			echo Detected: $ff is older than $script
			if grep $point $queue | grep -q $num; then 
				echo $nameExec for point $point at $num/$job is still running. Just wait.
			elif [ "$rep" == "true" ] ; then
				repair=(${repair[@]} "$num")
				echo Point $point will be repaired at $num/$job
			else
				echo $nameExec failed for point $point at $num/$job. You can resubmit with -f to repair
			fi
			continue
		fi

		# at this stage output file is newer than script
		if [ -s $log ] && ! tail $log | grep -q Finished ; then
			echo Detected: $log did not finish
			if grep $point $queue | grep -q $num; then 
				echo $nameExec for point $point at $num/$job is still running. Just wait.
			elif [ "$rep" == "true" ] ; then
				repair=(${repair[@]} "$num")
				echo Point $point will be repaired at $num/$job
			else
				echo $nameExec failed for point $point at $num/$job. You can resubmit with -f to repair
			fi
			continue
		fi
	done

	for rr in "${repair[@]}" ; do

		echo Repairing $point at $rr/$job

		this=$outdir/this_sensitivity.card
		if ! [ -s $this ] ; then
			# take cards from upper folder
			card=$outdir/../multi.card
			fitc=$outdir/../fit_options.card
			oscc=$outdir/../oscillation.card
			atmo=$outdir/../beam_sample.card
			beam=$outdir/../atmo_sample.card
			if ! [ -s $card ] ; then
				echo ERROR There is no main card $(realpath $card), very bad!
				exit 1
			fi

			sed -i "s:^point.*:point\t$t:" $card
			sed -i "s:^output.*:output\t\"$outdir/SpaghettiSens.root\":" $card

			cp $card $this
			sed -i "s:^fit_parameters.*:fit_parameters\t\"$fitc\":" $this
			sed -i "s:^oscillation_parameters.*:oscillation_parameters\t\"$oscc\":" $this
			sed -i "s:^beam_parameters.*:beam_parameters\t\"$beam\":" $this
			sed -i "s:^atmo_parameters.*:atmo_parameters\t\"$atmo\":" $this
		fi

		recover=$outdir/Recover.$rr.sub
		log=$outdir/L$nameExec.$num.log
		if [ "$SCHED" == "HTCONDOR" ] ; then
			cat > $recover << EOF
#! /bin/bash
# script submission for SLURM
# sumbit with --
#	$sub $recover

executable		= $Sens
arguments		= $rr $job $this
getenv			= True
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $log
error			= $log

queue
EOF
		elif [ "$SCHED" == "SLURM" ] ; then
			cat > $record << EOF
#! /bin/bash
# script submission for Slurm
# sumbit with --
#	$sub $recover

#SBATCH --job-name=fixing
#SBATCH -o $log
#SBATCH -p nms_research,shared
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1

srun $Sens $rr $job $this

EOF
		fi
		$sub $scriptname
	done
done

if [ "${#repeat[@]}" -gt 0 ] ; then
	point_file=".points_list"
	echo "${repeat[@]}" > .point_to_repeat
	card=$root/multi.card
	if ! [ -s $card ] ; then
		echo ERROR There is no main card $(realpath $card), very bad!
		exit 1
	fi

	$PWD/$trisens -x -r $root -p $point_file
fi
