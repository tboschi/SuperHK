#! /bin/bash

#oscillate files from GlobalOsc

Oscillator=/home/tboschi/OscAna/SuperHK/bin/oscillator
cardS=/home/tboschi/OscAna/SuperHK/cards
JM=/home/tboschi/jobManager
Queue=$JM/submitToQueue.sh
MAX_JOBS=500

global=/home/tboschi/data/
root=/home/tboschi/data/
MH=NH

#global contains already asim or hkdr

while getopts 'ins:g:' flag; do
	case "${flag}" in
		i) MH=IH ;;
		n) MH=NH ;;
		g) global="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

#where unoscillated files are
input=${global%/}/$MH/global
output=${global%/}/$MH/oscillated
nfiles=$(ls $input/*.root | wc -l)

echo $input
echo $output

#create files	#needs to be done once only

mkdir -p $output


cp $cardS/t2hkk_.card $cardS/t2hkk_oscillated.card

sed -i "s:input.*:input \"$input/*.root\":"	"$cardS/t2hkk_oscillated.card"
sed -i "s:output.*:output \"$output\":"		"$cardS/t2hkk_oscillated.card"

if [ $MH = "NH" ]; then
	sed -i "s:hierarchy.*:hierarchy \"normal\":"	"$cardS/t2hkk_oscillated.card"
else
	sed -i "s:hierarchy.*:hierarchy \"inverted\":"	"$cardS/t2hkk_oscillated.card"
fi

#oscillate beam part, updating cards with new systematics
count=0
chunk=$(expr $nfiles / 20)
rest=$(expr $nfiles % 20)
rm -f /home/tboschi/jobManager/*.list
rm -f $output/Roscillator_*.sh
for j in {0..19} ; do
	cp "$cardS/t2hkk_oscillated.card" $cardS/t2hkk_oscillated_$j".card" 
	first=$(expr $chunk \* $j)
	last=$(expr $first + $chunk - 1)
	echo $first $last 
	sed -i "s:first.*:first $first:" $cardS/t2hkk_oscillated_$j".card"
	sed -i "s:last.*:last $last:" $cardS/t2hkk_oscillated_$j".card"

	cat > $output/Roscillator_$j".sh" <<EOF
#!/bin/bash
source /home/tboschi/OscAna/Setup_OscAna.sh
$Oscillator $cardS/t2hkk_oscillated_$j.card
EOF
	echo -r Oscillator$j -o $output/Roscillator_$j.log -e $output/Roscillator_$j.err $output/Roscillator_$j".sh" >> /home/tboschi/jobManager/atmpd.list
done

if [ $rest -ne 0 ] ; then
	cp "$cardS/t2hkk_oscillated.card" $cardS/t2hkk_oscillated_20.card 
	first=$(expr $last + 1)
	last=$(expr $first + $rest)
	echo $first $last 
	sed -i "s:first.*:first $first:" $cardS/t2hkk_oscillated_20.card
	sed -i "s:last.*:last $last:" $cardS/t2hkk_oscillated_20.card
fi

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
	count=$(expr $count % ${#queues[@]})
	echo "count " $count
	mv $list $JM/${queues[$count]}.list
	list=$JM/${queues[$count]}.list
done

#mv $cardS/t2hkk_oscillated*.card $output/
