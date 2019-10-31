#! /bin/bash

#run GlobalOsc efficiently

Osc=/home/tboschi/OscAna/hk.atm+beam/submitOsc.sh
JM=/home/tboschi/jobManager
Queue=$JM/submitToQueue.sh
MAX_JOBS=500

global=""

#global already contains hkdr or asim
while getopts 'inf:g:' flag; do
	case "${flag}" in
		i) MH=IH ;;
		n) MH=NH ;;
		f) NFILES="${OPTARG}" ;;
		g) global="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

global=${global%/}

mkdir -p $global/$MH/global

cp $global/*.card $global/$MH/global/global.card
card=$global/$MH/global/global.card

#name=$(grep OutputOscFile $card | cut -d'"' -f2)
#output=${name%.*}

#echo sed -i "s:OutputOscFile.*:OutputOscFile \"$output.root\":" $card

if [[ $MH = "NH" ]]; then
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" $card
else
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" $card
fi

rm $JM/*.list

$Osc $card $global/$MH/global $NFILES

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
	count=$(expr $count % ${#queues[@]})
	echo "count " $count
	mv $list $JM/${queues[$count]}.list
	list=$JM/${queues[$count]}.list
done

echo "global finished"
