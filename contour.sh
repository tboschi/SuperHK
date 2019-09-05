#! /bin/bash

#Osc=/home/tboschi/OscAna/hk.atm+beam/submitOsc.sh

Sens=/home/tboschi/OscAna/hk.atm+beam/submitSens.sh
JM=/home/tboschi/jobManager
Queue=$JM/submitToQueue.sh
cardS=/home/tboschi/OscAna/SuperHK/cards
Exclude=/home/tboschi/OscAna/SuperHK/bin/exclusion
Penalty=/home/tboschi/OscAna/SuperHK/bin/addpenalty
Contour=/home/tboschi/OscAna/Osc3++/processing/build.contours/BuildContourPlots
#samples=(SK T2HK SKT2HK)
#samples=(T2HK SKT2HK)
samples=(T2HK)
MAX_JOBS=500

root=/home/tboschi/data/
global=""
MH_1=""
MH_2=""

while getopts '1:2:r:g:' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		2) MH_2="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
		g) global="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

global=${global%/}
root=${root%/}

NFILES=$(ls $global/$MH_1/oscillated/*.root | wc -l)

point=$(cat $global/point.info)
point=(${point})
mhfit=$MH_1"_"$MH_2


#penalty section

base=${root##*/}
cp $cardS/penalty_$base.card $root/$mhfit/sensitivity/penalty_sensitivity.card

for t in "${point[@]}" ; do
	rm -f $root/$mhfit/sensitivity/sens_t$t/*penalised*
done
sed -i "s:files.*:files \"$root/$mhfit/sensitivity/sens_t*/SpaghettiSens.*.root\":" $root/$mhfit/sensitivity/penalty_sensitivity.card
$Penalty $root/$mhfit/sensitivity/penalty_sensitivity.card
#rm -f $root/$mhfit/sensitivity/sens_t*/SpaghettiSens.*.??????.??????.root

for t in "${point[@]}" ; do
	mkdir -p $root/$mhfit/contours/sens_t$t
done

for ff in "${samples[@]}" ; do
	echo Processing $ff set
	for t in "${point[@]}" ; do
		#contour of pure files
		hadd $root/$mhfit/contours/sens_t$t/all.$ff.root $root/$mhfit/sensitivity/sens_t$t/SpaghettiSens.$ff.*.root
		$Contour $root/$mhfit/contours/sens_t$t/all.$ff.root >& /dev/null
		mv ChiSquared.root $root/$mhfit/contours/sens_t$t/uniform.$ff.root
		rm $root/$mhfit/contours/sens_t$t/all.$ff.root
		#contour of penalised files
		hadd $root/$mhfit/contours/sens_t$t/all.$ff.root $root/$mhfit/sensitivity/sens_t$t/SpaghettiSens.$ff'_penalised'.*.root
		$Contour $root/$mhfit/contours/sens_t$t/all.$ff.root >& /dev/null
		mv ChiSquared.root $root/$mhfit/contours/sens_t$t/gaussian.$ff.root
		rm $root/$mhfit/contours/sens_t$t/all.$ff.root
	done
done

if [ ${#point[@]} -gt 1 ]
then
	mkdir -p $root/$mhfit/exclusion/

	for ff in "${samples[@]}" ; do
		ls $root/$mhfit/sensitivity/sens_t*/SpaghettiSens.$ff'_penalised'.*.root > listExclusion
		$Exclude listExclusion	$root/$mhfit/exclusion/excl_$ff.dat
	done
fi

echo "DONE"
