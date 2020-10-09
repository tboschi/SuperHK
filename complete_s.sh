#! /bin/bash

Sens=$PWD/cross-fitter.sh
sub=sbatch

name=${1%.*}

while read -r point ; do
	# if finished is not shown, means something bad happened
	echo Checking point $point

	outdir=$name'_'$point

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
				cat > $scriptname << EOF
#! /bin/bash
# script submission for SLURM
# sumbit with --
#	$sub $scriptname

#SBATCH --job-name=fixing
#SBATCH -o $log
#SBATCH -p nms_research,shared
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1

srun $Sens $proc $all $outdir/this_sensitivity.card

EOF
				$sub $scriptname
			fi
		done
	fi
done < $1
