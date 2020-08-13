#! /bin/bash

pena=./penalise.sh
cont=./contour.sh
excl=./excludeall.sh

study="/data/tboschi/HKsens/OscAna/SuperHK/errorstudy"
point="asim/NH_NH"

#model=( "nuenorm1_corr" "nuenorm2_corr" "nuenorm3_corr" "nuenorm4_corr"  "nuenorm5_corr" 
	#"nuenorm1_anti" "nuenorm2_anti" "nuenorm3_anti" "nuenorm4_anti"  "nuenorm5_anti" )
#model=( "stats_untuned" )

#model=( "new0_25/asim/NH_NH" "new0_50/asim/NH_NH" "new0_75/asim/NH_NH" )
#model=( "new0_flux/asim/NH_NH" "flux_lukas_flux/asim/NH_NH" 
        #"12_/asim/NH_NH" "12a/asim/NH_NH" "12b/asim/NH_NH" )
#model=( "11a/asim/NH_NH" "11b/asim/NH_NH")
#model=( "nuenorm1_corr/asim/NH_NH" "nuenorm2_corr/asim/NH_NH" "nuenorm3_corr/asim/NH_NH" "nuenorm4_corr/asim/NH_NH" "nuenorm5_corr/asim/NH_NH" "nuenorm1_anti/asim/NH_NH" "nuenorm2_anti/asim/NH_NH" "nuenorm3_anti/asim/NH_NH" "nuenorm4_anti/asim/NH_NH" "nuenorm5_anti/asim/NH_NH" )
#model=( "stats/asim/NH_NH" "new_nuenorm5_corr/asim/NH_NH" "new_nuenorm5_anti/asim/NH_NH" )
#model=( "stats/asim/NH_NH" "new_nuenorm5_corr/asim/NH_NH" "0/asim/NH_NH" )
model=( "stats_fine/asim/NH_NH" "nominal_fine/asim/NH_NH" "nominal_energyscale_fine_2/asim/NH_NH" )

ff=false
while getopts 'f' flag; do
	case "${flag}" in
		f) ff=true ;;
		*) exit 1 ;;
	esac
done

for mod in "${model[@]}"
do
	#echo Penalise

	##$pena -r $study/$mod/$point/
	#$pena -r $study/$mod/

	#running=$(condor_q -all -format "%s\n" cmd | grep addpenalty | wc -l)
	#while [ $running -gt 0 ] ; do
	#	echo 'waiting 10s..'
	#	sleep 10
	#	running=$(condor_q -all -format "%s\n" cmd | grep addpenalty | wc -l)
	#done

	point=point_$(cat $study/$mod/sensitivity/point.info)
	echo Contours $point
	if [ "$ff" = true ]; then
		$excl -p -r $study/$mod/
	else
		$cont -p -r $study/$mod/ $point
	fi
done
