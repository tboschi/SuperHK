#! /bin/bash

sens=./trisens_c.sh

study="/data/tboschi/HKsens/OscAna/SuperHK/errorstudy"


#model=( "-s -r $study/stats/asim/ -1 NH -2 NH"
#	"-r $study/new_nuenorm5_corr/asim/ -1 NH -2 NH"
#model=(	"-r $study/0/asim/ -1 NH -2 NH")
model=(	"-r $study/test/ -1 NH -2 NH")

for mod in "${model[@]}" ; do
#while [ $count -lt $lengt ] ; do

	echo 'Analysing model' $mod

	echo $sens "$mod"
	$sens "$mod"

	while [ $(condor_q -run -format "%s\n" cmd | grep fitter | wc -l) -gt 50 ] ; do
		echo 'launcher: waiting for jobs to end...(2 min)'
		sleep 120
	done
done
