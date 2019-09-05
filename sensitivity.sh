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

mkdir -p $root/$mhfit/sensitivity/
cp $global/*.card $root/$mhfit/sensitivity/sensitivity.card
card=$root/$mhfit/sensitivity/sensitivity.card

#get ready to sensitivity

T2KnumuFile=$root/../systematics/FHC1Rmu.fij.t2k_p1.root
T2KnueFile=$root/../systematics/FHC1Re.fij.t2k_p1.root
T2KnumubarFile=$root/../systematics/RHC1Rmu.fij.t2k_p1.root
T2KnuebarFile=$root/../systematics/RHC1Re.fij.t2k_p1.root
matrix=$root/../systematics/combinedmatrix.root

sed -i "s:T2KnumuFijTable.*:T2KnumuFijTable \"$T2KnumuFile\":" $card
sed -i "s:T2KnueFijTable.*:T2KnueFijTable \"$T2KnueFile\":" $card
sed -i "s:T2KnumubarFijTable.*:T2KnumubarFijTable \"$T2KnumubarFile\":" $card
sed -i "s:T2KnuebarFijTable.*:T2KnuebarFijTable \"$T2KnuebarFile\":" $card
sed -i "s:AllErrFile.*:AllErrFile \"$matrix\":" $card

name=$(grep OutputOscFile $card | cut -d'"' -f2)
input=${name%.*}

echo $input 

card_1=${card/.card/_$MH_1"_1.card"}
card_2=${card/.card/_$MH_2"_2.card"}

cp $card $card_1
mv $card $card_2

sed -i "s:InputToFit.*:InputToFit \"$global/$MH_1/oscillated/$input.*.root\":" $card_1
sed -i "s:InputToFit.*:InputToFit \"$global/$MH_2/oscillated/$input.*.root\":" $card_2

if [[ $MH_1 = "NH" ]]; then
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" $card_1
else
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" $card_1
fi

if [[ $MH_2 = "NH" ]]; then
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" $card_2
else
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" $card_2
fi

#load fit
rm -f /home/tboschi/jobManager/*.list
for t in "${point[@]}" ; do
	$Sens $NFILES $card_1 $card_2 $t $root/$mhfit/sensitivity/sens_t$t
done

#smart submit to the queue system
queues=(atmpd ALL all lowe calib)
count=0
list=$JM/${queues[$count]}.list
echo 'Submitting to the queue system'
while [ -s $list ] ; do
	wcleft=$(wc -l $list)
	jobsinqueue=$(qstat ${queues[$count]} | grep run | awk '{print $1;}')
	echo "   jobs left to sumbit : " ${wcleft%%/*}
	echo "   jobs in " ${queues[$count]} " : " $jobsinqueue
	$Queue $(expr $MAX_JOBS - $jobsinqueue)
        count=$(expr $count + 1)

	#queues finished, put something on queue
	if [ $count -gt ${#queues[@]} ]
	then
		MAX_JOBS=$(expr $MAX_JOBS + 100)
	fi

	count=$(expr $count % ${#queues[@]})
	echo "count " $count
	mv $list $JM/${queues[$count]}.list
	list=$JM/${queues[$count]}.list
done

#wait until jobs are finished
#while true ; do
#	check=$(qstat -u tboschi | grep tboschi)
#	if [ -n "$check" ]; then
#		echo 'waiting 5min...'
#		sleep 300
#	else
#		break
#	fi
#done


#penalty section

#base=${root##*/}
#cp $cardS/penalty_$base.card $root/$mhfit/sensitivity/penalty_sensitivity.card
#
#sed -i "s:files.*:files \"$root/$mhfit/sensitivity/sens_t*/SpaghettiSens.*.root\":" $root/$mhfit/sensitivity/penalty_sensitivity.card
#$Penalty $root/$mhfit/sensitivity/penalty_sensitivity.card
##rm -f $root/$mhfit/sensitivity/sens_t*/SpaghettiSens.*.??????.??????.root
#
#for t in "${point[@]}" ; do
#	mkdir -p $root/$mhfit/contours/sens_t$t
#done
#
#for ff in "${samples[@]}" ; do
#	echo Processing $ff set
#	for t in "${point[@]}" ; do
#		hadd $root/$mhfit/contours/sens_t$t/all.$ff.root $root/$mhfit/sensitivity/sens_t$t/SpaghettiSens.$ff'_penalised'.*.root
#		$Contour $root/$mhfit/contours/sens_t$t/all.$ff.root >& /dev/null
#		mv ChiSquared.root $root/$mhfit/contours/sens_t$t/contour.$ff.root
#		rm $root/$mhfit/contours/sens_t$t/all.$ff.root
#	done
#done
#
#if [ ${#point[@]} -gt 1 ]
#then
#	mkdir -p $root/$mhfit/exclusion/
#
#	for ff in "${samples[@]}" ; do
#		ls $root/$mhfit/sensitivity/sens_t*/SpaghettiSens.$ff'_penalised'.*.root > listExclusion
#		$Exclude listExclusion	$root/$mhfit/exclusion/excl_$ff.dat
#	done
#fi
#
#echo "DONE"
