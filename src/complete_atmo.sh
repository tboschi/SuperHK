#! /bin/bash

usage="usage: $0 [<options>] <dir>

Check if jobs were completed and repair if needed where <dir> is
the output folder with the precompute atmospheric input files.
Works with HTCondor, Slurm, SGE or PBS.

Options
    -f	      repair files, otherwise just display information
    -i	      if triatmo is running, just ignore it
    -w <dir>  if the log files are in a different directory,
              please specify it
    -h        print this message
"

rep=false
skip=false
logr=""
while getopts 'i:r:w:fh' flag; do
	case "${flag}" in
		f) rep=true ;;
		i) skip=true ;;
		w) logr="${OPTARG}" ;;
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
	echo -e "\n    $0 -h\n"
	exit 1
elif [ ! -d "$1" ] ; then
	echo Argument \"$1\" is not valid. Check usage with
	echo -e "\n    $0 -h\n"
	exit 1
fi

Bin=$PWD/bin/cross-binary
Atmo=atmo_input
triatmo=bin/triatmo


if [ "$skip" == "false" ] && pgrep $triatmo ; then
	echo $triatmo is still running. Wait for it to finish.
	exit 1
fi


SCHED=""
queue=".queue_status"
if condor_q &> /dev/null ; then
	SCHED="HTCONDOR"
	sub=condor_submit
	generate=src/slurm.recover
	condor_q -autoformat cmd args > $queue
elif squeue &> /dev/null ; then
	SCHED="SLURM"
	sub=sbatch
	generate=src/slurm.recover
	squeue -h -r -u $USER -o "%o %K" > $queue
elif qstat -Q &> /dev/null ; then
	SCHED="PBS"
	sub=qsub
	generate=src/sge.recover
	qstat -u $USER -f > $queue
elif qstat -sql &> /dev/null ; then
	SCHED="SGE"
	sub=qsub
	generate=src/sge.recover
	qstat -u $USER -f > $queue
else
	echo "There is neither HTCondor nor Slurm nor qsub \(SGE or PBS\) on this machine. \
		Iam sorry, I cannot help you"
	exit 1
fi

outdir=$(realpath $1)
script=$(ls -rt $outdir/Ratmo_input.sub | tail -n1)

if [ -n "$logr" ] ; then
	outlog=$logr
else
	outlog=$outdir
fi

# find all root files
all=$(find $outdir -name "atmo.*.root")
all=(${all})

# find how many jobs there should be
if grep -q queue $script ; then
	njobs=$(grep queue $script | cut -f2 -d' ')
elif grep -q srun $script ; then
	njobs=$(grep srun $script | cut -f5 -d' ')
elif grep -q "PBS" $script ; then
	njobs=$(tail -n1 $script | cut -f4 -d' ')
else	# it is SGE
	njobs=$(tail -n1 $script | cut -f4 -d' ')
fi
echo "expected jobs " $njobs

# now all known files
for out in "${all[@]}" ; do

	# skip file with weird format
	if ! [[ $out =~ atmo\.[0-9]+\.root ]] ; then
		echo Detected: skip unknown file $out
		continue
	fi

	num=${out%.root}
	num=${num##*.}

	if [ $num -ge $njobs ] ; then
		echo Detected: file misplaced $out 
		if [ "$rep" == "true" ] ; then
			rm $out
		else
			echo "    " you can resubmit with -f to delete
		fi
		continue
	fi
done

repair=()
for num in $(seq 0 $((njobs - 1)) ) ; do 

	out=$outdir/atmo.$num.root
	log=$outlog/Latmo_input.$num.log

	bad=false
	if ! [ -s $out ] ; then	
		echo Detected: $out does not exist
		bad=true
	elif [ $out -ot $script ] ; then
		echo Detected: $out is older than $script
		bad=true  # at this stage output file is newer than script
	elif ! [ -s $log ] ; then
		echo Detected: $log does not exist
		bad=true
	elif  ! tail -n1 $log | grep -q Finished ; then
		echo Detected: $log did not finish. Last lines are
		tail "$log"
		echo ""
		bad=true
	fi
	#elif head $log | grep -q permission ; then
	#	echo Detected: $log did not finish. Last lines are
	#	bad=true

	if [ "$bad" == "true" ] ; then
		if grep -q $num $queue; then 
			echo $Atmo for point $point at $num/$njobs is still running. Just wait.
		elif [ "$rep" == "true" ] ; then
			repair=(${repair[@]} "$num")
			echo Point $point will be repaired at $num/$njobs
		else
			echo $Atmo failed for point $point at $num/$job. You can resubmit with -f to repair
		fi
	fi
done

this=$outdir/multi.card
if [ "$rep" == "true" ] && ! [ -s $this ] ; then
	echo "Missing card file, cannot repair!"
	echo "Please repair manually or resubmit full $triatmo"
	exit 1
fi

for rr in "${repair[@]}" ; do

	echo Repairing $rr/$njobs

	recover=$outdir/Recover.$rr.sub

	$generate $Bin $Atmo $njobs $this $outlog $rr > $recover
	$sub $recover
done

#	if [ "$SCHED" == "HTCONDOR" ] ; then
#		cat > $recover << EOF
##! /bin/bash
## script submission for condor
## sumbit with --
##	$sub $recover
#
#executable		= $Atmo
#arguments		= $rr $job $this
#getenv			= True
#should_transfer_files	= IF_NEEDED
#when_to_transfer_output	= ON_EXIT
#initialdir		= $PWD
#output			= $outlog/L$nameExec.$rr.log
#error			= $outlog/L$nameExec.$rr.log
#stream_output		= True
#stream_error		= True
#+JobFlavor="testmatch"
#
#queue
#
#EOF
#		elif [ "$SCHED" == "SLURM" ] ; then
#			cat > $recover << EOF
##! /bin/bash
## script submission for Slurm
## sumbit with --
##	$sub $recover
#
##SBATCH --job-name=fixing
##SBATCH -o $log
##SBATCH -p nms_research,shared
##SBATCH --time=3-0
##SBATCH --cpus-per-task=1
#
#srun $Atmo fitter $rr $job $this
#
#EOF
#		fi
