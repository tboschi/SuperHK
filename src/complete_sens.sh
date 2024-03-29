#! /bin/bash

usage="usage: $0 [<options>] <list>

Check if jobs were completed and repair if needed where <list> is
the file with the list of fitted points in the sensitivity folder.
Works with HTCondor, Slurm, SGE or PBS.

Options
    -f	      repair files, otherwise just display information
    -i	      if trisens is running, just ignore it
    -w <dir>  if the log files are in a different directory,
              please specify it
    -h        print this message
"

rep=false
skip=false
logr=""
while getopts 'ifhw:' flag; do
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
elif [ ! -s "$1" ] ; then
	echo Argument \"$1\" is not valid. Check usage with
	echo -e "\n    $0 -h\n"
	exit 1
fi

Bin=$PWD/bin/cross-binary
Sens=fitter
trisens=bin/trisens


if [ "$skip" == "false" ] && pgrep $trisens ; then
	echo $trisens is still running. Wait for it to finish.
	exit 1
fi


SCHED=""
queue=".queue_status"
if condor_q &> /dev/null ; then
	SCHED="HTCONDOR"
	sub=condor_submit
	generate=src/htcondor.recover
	condor_q -autoformat cmd args > $queue
elif squeue &> /dev/null ; then
	SCHED="SLURM"
	sub=sbatch
	generate=src/slurm.recover
	squeue -h -r -u $USER -o "%o %K" > $queue
elif qstat &> /dev/null ; then
	SCHED="SGE"
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

name=${1%.*}
root=${name%/*}

repeat=()
while read -r point ; do
	# if finished is not shown, means something bad happened
	echo Checking point $point

	outdir=$name'_'$point
	outdir=$(realpath $outdir)
	script=$(ls -rt $outdir/R*.$point.sh | tail -n1)

	if [ -n "$logr" ] ; then
		outlog=$logr/$(expr "$outdir" : '.*\(.H_.H/.*\)')
	else
		outlog=$outdir
	fi

	# directory does not exists or no script
	if ! [ -d "$outdir" ] || ! [ -s "$script" ] ; then
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
	all=$(find $outdir -name "SpaghettiSens.*.root")
	all=(${all})

	if [ "${#all[@]}" -eq 0 ] ; then # no files
		# check if job is still running
		echo Detected: there are no output files in $outdir
		if grep -q $point $queue; then 
			echo $Sens for point $point is still running. Just wait.
		elif [ "$rep" == "true" ] ; then
			echo Point $point will be resubmitted
			repeat=(${repeat[@]} $point)
		else
			echo $Sens for point $point failed and never ran. You can resubmit with -f to repair
		fi
		continue
	fi 


	# find how many jobs there should be
	if grep -q queue $script ; then
		njobs=$(grep queue $script | cut -f2 -d' ')
	elif grep -q srun $script ; then
		njobs=$(grep srun $script | cut -f5 -d' ')
	elif grep -q "PBS" $script ; then
		njobs=$(tail -n1 $script | cut -f4 -d' ')
	else	# it is SGE
		njobs=$(tail -n2 $script | cut -f4 -d' ')
	fi

	# now all known files
	for out in "${all[@]}" ; do

		# skip file with weird format
		if ! [[ "$out" =~ SpaghettiSens\.[0-9]+\.root ]] ; then
			echo Detected: skip unknown file $out
			continue
		fi

		num=${out%.root}
		num=${num##*.}
		
		if [ "$num" -ge "$njobs" ] ; then
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

		out=$outdir/SpaghettiSens.$num.root
		log=$outlog/L$Sens.$num.log

		bad=false
		if ! [ -s "$out" ] ; then	
			echo Detected: $out does not exist
			bad=true
		elif [ "$out" -ot "$script" ] ; then
			echo Detected: $out is older than $script
			bad=true  # at this stage output file is newer than script
		elif ! [ -s "$log" ] ; then
			echo Detected: $log does not exist
			bad=true
		elif ! tail -n1 "$log" | grep -q Finished ; then
			echo Detected: $log did not finish. Last lines are
			tail "$log"
			echo ""
			bad=true
		fi

		if [ "$bad" == "true" ] ; then
			if grep $point $queue | grep -q $num; then 
				echo $Sens for point $point at $num/$njobs is still running. Just wait.
			elif [ "$rep" == "true" ] ; then
				repair=(${repair[@]} "$num")
				echo Point $point will be repaired at $num/$njobs
			else
				echo $Sens failed for point $point at $num/$job. You can resubmit with -f to repair
			fi
		fi
	done

	for rr in "${repair[@]}" ; do

		echo Repairing $point at $rr/$njobs

		this=$outdir/this_sensitivity.card
		if ! [ -s "$this" ] ; then
			# take cards from upper folder
			card=$outdir/../multi.card
			fitc=$outdir/../fit_options.card
			oscc=$outdir/../oscillation.card
			atmo=$outdir/../beam_sample.card
			beam=$outdir/../atmo_sample.card
			if ! [ -s "$card" ] ; then
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


		recover=$outdir/Recover.$rr.sh
		#log=$outlog/L$nameExec.$rr.log

		$generate $Bin $Sens $njobs $this $outlog $rr > $recover
		$sub $recover
	done
done < $1

if [ "${#repeat[@]}" -gt 0 ] ; then
	point_file=".points_list"
	echo "${repeat[@]}" > $point_file
	card=$root/multi.card
	if ! [ -s "$card" ] ; then
		echo ERROR There is no main card $(realpath $card), very bad!
		exit 1
	fi

	scan=""
	if [[ "$1" =~ .*CPV.* ]] ; then
		scan="-f CPV"
	fi

	$PWD/$trisens -x $scan -r $root -p $point_file
fi


#		if [ "$SCHED" == "HTCONDOR" ] ; then
#			cat > $recover << EOF
##! /bin/bash
## script submission for HTCondor
## sumbit with --
##	$sub $recover
#
#executable		= $Sens
#arguments		= fitter $rr $job $this
#getenv			= True
#should_transfer_files	= IF_NEEDED
#when_to_transfer_output	= ON_EXIT
#initialdir		= $PWD
#output			= $log
#error			= $log
#
#queue
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
##SBATCH --time=4-0
##SBATCH --cpus-per-task=1
##SBATCH --exclude=nodek50,nodek33
#
#srun $Sens fitter $rr $job $this
#
#EOF
#		fi
