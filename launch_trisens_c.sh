#! /bin/bash

sens=./trisens_c.sh

study="/data/tboschi/HKsens/OscAna/SuperHK/errorstudy"
point="asim"

#model=("0" "11a" "11b" "8")

#model=( "nuenorm1_corr" "nuenorm2_corr" "nuenorm3_corr" "nuenorm4_corr"  "nuenorm5_corr" 
	#"nuenorm1_anti" "nuenorm2_anti" "nuenorm3_anti" "nuenorm4_anti"  "nuenorm5_anti" )

model=( "nuenorm1_corr" "nuenorm2_corr" "nuenorm3_corr" "nuenorm4_corr" "nuenorm5_corr"
        "nuenorm1_anti" "nuenorm2_anti" "nuenorm3_anti" "nuenorm4_anti" "nuenorm5_anti" )

ff=false
mm=""
while getopts 'f' flag; do
	case "${flag}" in
		f) ff=true ;;
		*) exit 1 ;;
	esac
done

if [ "$ff" = true ]; then
	fop="-f"
else
	fop=""
fi

#shift $((OPTIND-1))

#while [ $(condor_q -run -format "%s\n" cmd | wc -l) -gt 50 ] ; do
#	echo 'launcher: waiting for jobs to end...(5 min)'
#	sleep 300
#done

for mod in "${model[@]}" ; do
#while [ $count -lt $lengt ] ; do

	echo 'Analysing model' $mod

	if [[ $mod = *"indy"* ]]; then
		mm="-m identity"
	else
		mm=""
	fi

	stats=${mod##*_}

	echo $sens $fop -g $study/global/$point -r $study/$mod/$point -1 NH -2 NH
	$sens $fop -g $study/global/$point -r $study/$mod/$point -1 NH -2 NH

	while [ $(condor_q -run -format "%s\n" cmd | grep fitter | wc -l) -gt 50 ] ; do
		echo 'launcher: waiting for jobs to end...(2 min)'
		sleep 120
	done
done
