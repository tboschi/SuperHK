#! /bin/bash

pena=./penalise.sh
cont=./contour.sh
excl=./excludeall.sh

study="/data/tboschi/HKsens/errorstudy"
point="asim/NH_NH"

#model=( "0" "11a" "11b" "8" )
#model=( "nuenorm1_corr" "nuenorm2_corr" "nuenorm3_corr" "nuenorm4_corr"  "nuenorm5_corr" 
	#"nuenorm1_anti" "nuenorm2_anti" "nuenorm3_anti" "nuenorm4_anti"  "nuenorm5_anti" )
#model=( "stats_untuned" )
model=( "0_E" "0_M" )

ff=false
while getopts 'f' flag; do
	case "${flag}" in
		f) ff=true ;;
		*) exit 1 ;;
	esac
done

for mod in "${model[@]}"
do
	echo Penalise
	$pena -r $study/$mod/$point/

	running=$(condor_q -all -submitter tboschi -format "%s\n" ProcId | wc -l)
	while [ $running -gt 0 ] ; do
		echo 'waiting 1min...'
		sleep 60
		running=$(condor_q -all -submitter tboschi -format "%s\n" ProcId | wc -l)
	done

	echo Contours
	if [ "$ff" = true ]; then
		$excl -p -r $study/$mod/$point/
	else
		$cont -p -r $study/$mod/$point/ point_49773
	fi
done
