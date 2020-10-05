#! /bin/bash

Sens=$PWD/cross-fitter.sh

root=$PWD/errorstudy
while getopts 'r:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done

# add PWD 
if [[ "$root" != /* ]]
then
	root=$PWD/${root%/}
else
	root=${root%/}
fi

find $root -name "*.log" > alllogs.list

while read -r log ; do
	# if finished is not shown, means something bad happened
	if ! tail $log | grep -q Finished ; then
		echo $log not ok
		#log is #/long/path/to/folder/CPV_scan_point/fitter.proc.log
		proc=${log##*/}
		#proc is fitter.proc.log
		proc=${proc#*.}
		#proc is proc.log
		proc=${proc%.*}
		#proc is proc

		outdir=${log%/*}
		#outidr is #/long/path/to/folder/CPV_scan_point/
		all=$(ls $outdir/*.log | wc -l)

		echo $outdir and $proc
		scriptname=$outdir/Recover.$proc.sub
		cat > $scriptname << EOF
# script submission for condor
# sumbit with --
#	condor_submit $scriptname

executable		= $Sens
arguments		= $proc $all $outdir/this_sensitivity.card
getenv			= True
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $log
error			= $log

queue
EOF
		echo condor_submit $scriptname
	fi
done < "alllogs.list"
