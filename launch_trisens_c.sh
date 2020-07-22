#! /bin/bash

sens=./trisens_c.sh

study="/data/tboschi/HKsens/OscAna/SuperHK/errorstudy"
point="asim"

model=("0") #"11a" "11b" "8")

#model=( "nuenorm1_corr" "nuenorm2_corr" "nuenorm3_corr" "nuenorm4_corr"  "nuenorm5_corr" 
	#"nuenorm1_anti" "nuenorm2_anti" "nuenorm3_anti" "nuenorm4_anti"  "nuenorm5_anti" )

#model=( "-s -r $study/stats/asim/ -1 NH -2 NH"
#	"-r $study/new_nuenorm5_corr/asim/ -1 NH -2 NH"
#model=(	"-r $study/0/asim/ -1 NH -2 NH")
model=(	"-r $study/nominal_energyscale_fine_2/asim/ -1 NH -2 NH")

for mod in "${model[@]}" ; do
#while [ $count -lt $lengt ] ; do

	echo 'Analysing model' $mod

	echo $sens -g $study/global/$point "$mod"
	$sens -g $study/global/$point $mod

	while [ $(condor_q -run -format "%s\n" cmd | grep fitter | wc -l) -gt 50 ] ; do
		echo 'launcher: waiting for jobs to end...(2 min)'
		sleep 120
	done
done
