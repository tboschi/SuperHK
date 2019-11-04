#! /bin/bash

#Osc=/home/tboschi/OscAna/hk.atm+beam/submitOsc.sh

Sens=/home/tboschi/OscAna/hk.atm+beam/submitSens.sh
Stat=/home/tboschi/OscAna/hk.atm+beam/submitSensStatOnly.sh
Fast=/home/tboschi/OscAna/hk.atm+beam/submitSensFast.sh
Fastat=/home/tboschi/OscAna/hk.atm+beam/submitSensStatFast.sh
JM=/home/tboschi/jobManager
Queue=$JM/simpleSubmit.sh
cardS=/home/tboschi/OscAna/SuperHK/cards
samples=(T2HK)

root=/home/tboschi/data/
global=""
MH_1=""
MH_2=""
ss=false	#false if systematic fit, true if only stat fit
ff=false	#false if systematic fit, true if fast fit for dCP scan
mtype="correlation"

while getopts '1:2:r:g:m:sf' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		2) MH_2="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
		g) global="${OPTARG}" ;;
		s) ss=true ;;
		f) ff=true ;;
		m) mtype="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

if [ "$ss" = true ] && [ "$ff" = true ]
then
	pinfo="scan.info"
	tname="scan_"
	Exec=$Fastat
elif [ "$ss" = true ] && [ "$ff" = false ]
then
	pinfo="point.info"
	tname="point_"
	Exec=$Stat
elif [ "$ss" = false ] && [ "$ff" = true ]
then
	pinfo="scan.info"
	tname="scan_"
	Exec=$Fast
else
	pinfo="point.info"
	tname="point_"
	Exec=$Sens
fi
	

global=${global%/}
root=${root%/}

NFILES=$(ls $global/$MH_1/oscillated/*.root | wc -l)

mhfit=$MH_1"_"$MH_2

mkdir -p $root/$mhfit/sensitivity/
cp $global/*.card $root/$mhfit/sensitivity/sensitivity.card
card=$root/$mhfit/sensitivity/sensitivity.card

if [ ! -s $global/$pinfo ]
then
	echo "Point file does not exist"
fi

point=$(cat $global/$pinfo)
point=(${point})

#get ready to sensitivity
#edit card to point to correct systematics
T2KnumuFile=$root/../systematics/FHC1Rmu.fij.t2k_spline.root
T2KnueFile=$root/../systematics/FHC1Re.fij.t2k_spline.root
T2KnumubarFile=$root/../systematics/RHC1Rmu.fij.t2k_spline.root
T2KnuebarFile=$root/../systematics/RHC1Re.fij.t2k_spline.root
matrix=$root/../systematics/combinedmatrix.root
#mtype="correlation"

sed -i "s:T2KnumuFijTable.*:T2KnumuFijTable \"$T2KnumuFile\":" $card
sed -i "s:T2KnueFijTable.*:T2KnueFijTable \"$T2KnueFile\":" $card
sed -i "s:T2KnumubarFijTable.*:T2KnumubarFijTable \"$T2KnumubarFile\":" $card
sed -i "s:T2KnuebarFijTable.*:T2KnuebarFijTable \"$T2KnuebarFile\":" $card
sed -i "s:AllErrFile.*:AllErrFile \"$matrix\":" $card
sed -i "s:AllErrCov.*:AllErrCov \"$mtype\":" $card

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

MAX_JOBS=300
MAX_QUEUE=100
queues=(atmpd calib ALL all lowe)
#load fit
for t in "${point[@]}" ; do
	#create new jobs list
	rm -f /home/tboschi/jobManager/*.list
	$Exec $NFILES $card_1 $card_2 $t $root/$mhfit/sensitivity/$tname$t

	#smart submit to the queue system
	count=0
	list=$JM/${queues[$count]}.list
	echo 'Submitting to the queue system'
	while [ -s $list ] ; do
		wcleft=$(wc -l $list)
		jobsrunning=$(qstat -a ${queues[$count]} | head -n2 | grep run | awk '{print $1;}')
		jobsinqueue=$(qstat -a ${queues[$count]} | head -n2 | grep run | awk '{print $3;}')
		jobstot=$(expr $jobsrunning + $jobsinqueue)
		echo "jobs left to sumbit : " ${wcleft%%/*}
		echo "jobs in " ${queues[$count]} " : " $jobsrunning " running and " $jobsinqueue " in queue"
		if [ $MAX_JOBS  -gt $jobsrunning ]
		then
			$Queue $(expr $MAX_JOBS  - $jobsrunning)
		fi
		if [ $MAX_QUEUE -gt $jobsinqueue ]
		then
			$Queue $(expr $MAX_QUEUE - $jobsinqueue)
		fi
		#if [ $(expr $MAX_JOBS - $jobsrunning) -gt 0 ]; then
		#	$Queue $(expr $MAX_JOBS - $jobsrunning)
		#fi
		count=$(expr $count + 1)
		echo "count " $count

		count=$(expr $count % ${#queues[@]})
		#queues finished, put something on queue
		#if [ $count -eq 0 ]
		#then
		#	MAX_JOBS=$(expr $MAX_JOBS + 100)
		#fi
		#echo "Max jobs " $MAX_JOBS

		mv $list $JM/${queues[$count]}.list
		list=$JM/${queues[$count]}.list
	done

	#wait until jobs are finished before moving to next point

	while [ $t -ne ${point[${#point[@]}-1]} ] ; do
		check=$(qstat -u tboschi | grep tboschi | wc -l)
		if [ $check -ge 300 ]; then
			echo 'waiting 10min...'
			sleep 600
		else
			break
		fi
	done
done
