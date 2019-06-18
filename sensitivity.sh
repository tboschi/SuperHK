#! /bin/bash

#Osc=/home/tboschi/OscAna/hk.atm+beam/submitOsc.sh

Sens=/home/tboschi/OscAna/hk.atm+beam/submitSens.sh
JM=/home/tboschi/jobManager
Queue=$JM/submitToQueue.sh
cardS=/home/tboschi/OscAna/SuperHK/cards
Exclude=/home/tboschi/OscAna/SuperHK/bin/exclusion
Penalty=/home/tboschi/OscAna/SuperHK/bin/addpenalty
Contour=/home/tboschi/OscAna/Osc3++/processing/build.contours/BuildContourPlots


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
mhfit=$MH_1"_"$MH_2

mkdir -p $root/$mhfit/sensitivity/
cp $global/*.card $root/$mhfit/sensitivity/sensitivity.card
card=$root/$mhfit/sensitivity/sensitivity.card

#get ready to sensitivity

T2KnumuFile=$root/../systematics/FHC1Rmu.fij.t2k_p1.root
T2KnueFile=$root/../systematics/FHC1Re.fij.t2k_p1.root
T2KnumubarFile=$root/../systematics/RHC1Rmu.fij.t2k_p1.root
T2KnuebarFile=$root/../systematics/RHC1Re.fij.t2k_p1.root

echo sed -i "s:T2KnumuFijTable.*:T2KnumuFijTable $T2KnumuFile:" $card
echo sed -i "s:T2KnueFijTable.*:T2KnueFijTable $T2KnueFile:" $card
echo sed -i "s:T2KnumubarFijTable.*:T2KnumubarFijTable $T2KnumubarFile:" $card
echo sed -i "s:T2KnuebarFijTable.*:T2KnuebarFijTable $T2KnuebarFile:" $card

name=$(grep OutputOscFile $card | cut -d'"' -f2)
input=${name%.*}

echo $input 

card_1=${card/.card/_$MH_1"_1.card"}
card_2=${card/.card/_$MH_2"_2.card"}

echo cp $card $card_1
echo mv $card $card_2

echo sed -i "s:InputToFit.*:InputToFit \"$global/$mhfit/oscillated/$input.*.root\":" $card_1
echo sed -i "s:InputToFit.*:InputToFit \"$global/$mhfit/oscillated/$input.*.root\":" $card_2

if [[ $MH_1 = "NH" ]]; then
	echo sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" $card_1
else
	echo sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" $card_1
fi

if [[ $MH_2 = "NH" ]]; then
	echo sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" $card_2
else
	echo sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" $card_2
fi

#perform fit
echo $Sens $NFILES $card_1 $card_2 $point $root/$mhfit/sensitivity/sens_t$point
#point=0
#for i in {0..100}
#do
#        point=$(expr 343 \* $i + 25)
#	echo test point $point
#	echo $Sens $NFILES $cardir/temp1.card $cardir/temp2.card $point $root/sensitivity/$mhfit/sens_t$point
#done

#smart submit to the queue system
queues=(atmpd ALL all lowe calib)
count=0
list=$JM/${queues[$count]}.list
echo 'Submitting to the queue system'
while [ -s $list ] ; do
	wcleft=$(wc -l $list)
	echo "   jobs left to sumbit : " ${wcleft%%/*}
	echo "   -> to queue " $list
	echo $Queue
        count=$(expr $count + 1)
	echo mv $list $JM/${queues[$count]}.list
	list=$JM/${queues[$count]}.list
done


#wait until jobs are finished
#while true ; do
while false ; do
	check=$(qstat -u tboschi | grep tboschi)
	if [[ $check ]]; then
		echo 'waiting 1min...'
		echo sleep 60
	else
		break
	fi
done

if   [[ $root == *"asim"* ]]; then
	cp $cardS/penalty_asimov.card $root/$mhfit/sensitivity/
elif [[ $root == *"hkdr"* ]]; then
	cp $cardS/penalty_hkdr $root/$mhfit/sensitivity/
fi

sed -i "s:files.*:files \"$root/$mhfit/sensitivity/sens_t*/SpaghettiSens.SK.*.root\":" $root/$mhfit/sensitivity/penalty_sensitivity.card
echo $Penalty $root/$mhfit/sensitivity/penalty_sensitivity.card
sed -i "s:files.*:files \"$root/$mhfit/sensitivity/sens_t*/SpaghettiSens.T2HK.*.root\":" $root/$mhfit/sensitivity/penalty_sensitivity.card
echo $Penalty $root/$mhfit/sensitivity/penalty_sensitivity.card
sed -i "s:files.*:files \"$root/$mhfit/sensitivity/sens_t*/SpaghettiSens.SKT2HK.*.root\":" $root/$mhfit/sensitivity/penalty_sensitivity.card
echo $Penalty $root/$mhfit/sensitivity/penalty_sensitivity.card

ls $root/$mhfit/sensitivity/sens_t*/SpaghettiSens_corrected.SK.*.root	  > list_SK
ls $root/$mhfit/sensitivity/sens_t*/SpaghettiSens_corrected.T2HK.*.root	  > list_T2HK
ls $root/$mhfit/sensitivity/sens_t*/SpaghettiSens_corrected.SKT2HK.*.root > list_SKT2HK

mkdir -p $root/$mhfit/contours/

files=(SK T2HK SKT2HK)
for ff in "${files[@]}" ; do
	echo Processing $ff set
	list="list_"$ff
	while IFS='' read -r file || [[ -n "$file" ]] ; do
		echo "Contouring: $file"
		echo "$Contour $file	&> /dev/null"
		file=${file/SpaghettiSens.$ff/contour.$ff}
		file=${file/.0*./.}
		echo saving to $file
		mv ChiSquared.root $file
	done < "$list"
	echo hadd $root/$mhfit/contours/contour.$ff'.'*'.root' $root/$mhfit/contours/contour.$ff'.root'
done

mkdir -p $root/$mhfit/exclusion/

echo $Exclude list_SK	  $root/$mhfit/exclusion/excl_SK.dat
echo $Exclude list_T2HK	  $root/$mhfit/exclusion/excl_T2HK.dat
echo $Exclude list_SKT2HK $root/$mhfit/exclusion/excl_SKT2HK.dat

echo "DONE"
