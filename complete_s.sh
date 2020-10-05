#! /bin/bash

Sens=$PWD/cross-fitter.sh
sub=sbatch

while getopts 'l:' flag; do
	case "${flag}" in
		l) list="${OPTARG}" ;;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done

name=${list%.*}

while read -r point ; do
	# if finished is not shown, means something bad happened
	echo checking point $point

	outdir=$name'_'$point

	allog=$(find $outdir -name "*.log")
	allog=(${allog})

	if [ "${#allog[@]}" -eq 0 ]; then
		echo No log file for point $point, job did not even start
		if [ -s "$outdir/R"*".$point.sub" ] ; then
			echo $sub $outdir/R*.$point.sub
		fi
	else
		for log in "${allog[@]}" ; do
			if ! tail $log | grep -q Finished ; then
				echo $point not ok
				#log is #/long/path/to/folder/CPV_scan_point/fitter.proc.log
				proc=${log##*/}
				#proc is fitter.proc.log
				proc=${proc#*.}
				#proc is proc.log
				proc=${proc%.*}
				#proc is proc

				all=${#allog[@]}

				echo $outdir and $proc
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

echo srun $Sens $proc $all $outdir/this_sensitivity.card

EOF
				echo $sub $scriptname
			fi
		done
	fi
done < $list
