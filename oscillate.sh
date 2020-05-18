#! /bin/bash

#oscillate files from GlobalOsc

Exec=/data/tboschi/HKsens/OscAna/SuperHK/bin/oscillator
card=/data/tboschi/HKsens/OscAna/SuperHK/cards

global=/home/tboschi/data/
root=/home/tboschi/data/

#global contains already asim or hkdr

while getopts 'g:' flag; do
	case "${flag}" in
		g) global="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

#where unoscillated files are
input=${global%/}/global
output=${global%/}/oscillated
#NFILES=$(ls $input/*.root | wc -l)
NFILES=307

echo $input
mkdir -p $output

cp $card/t2hkk_.card $card/t2hkk_oscillated.card

sed -i 's:input.*:input \"$input/*.root\":'	$card'/t2hkk_oscillated.card'
sed -i 's:output.*:output \"$output\":'		$card'/t2hkk_oscillated.card'

ih=$(grep InvertedHierarchy $input/global.card | cut -d ' ' -f 2)

if [ $ih -eq '0' ]; then
	sed -i 's:hierarchy.*:hierarchy \"normal\":'	$card/'t2hkk_oscillated.card'
elif [ $ih -eq '1' ]; then                                     
	sed -i 's:hierarchy.*:hierarchy \"inverted\":'	$card/'t2hkk_oscillated.card'
fi

nameExec=${Exec##*bin/}
rm -f $output/R$nameExec*.sh

#oscillate beam part, updating cards with new systematics
N_JOBS=20
chunk=$(expr $NFILES / $N_JOBS)
resto=$(expr $NFILES % $N_JOBS)

echo "ciao"
echo $chunk, $resto

if [ $resto -ne '0' ] ; then
	N_JOBS=$(expr $N_JOBS + 1)
fi

scriptname=$output/'R'$nameExec'.sh'
MAX_JOBS=100

cat > $scriptname << EOF
#$ -S /bin/bash
#$ -N L$nameExec
#$ -t 1-$N_JOBS
#$ -tc $MAX_JOBS
#$ -o $output/\$JOB_NAME.\$TASK_ID.log
#$ -j y

fnum=\$(expr \$SGE_TASK_ID - 1)

first=\$(expr $chunk \* \$fnum)
if [ \$SGE_TASK_ID -eq \$SGE_TASK_LAST ] ; then
	last=\$(expr \$first + $resto)
else
	last=\$(expr \$first + $chunk - 1)
fi	

echo Run \$fnum : oscillating files from set \$first to set \$last 
sed -i 's:first.*:first \$first:' '$card/t2hkk_oscillated_'\$fnum'.card'
sed -i 's:last.*:last \$last:'	  '$card/t2hkk_oscillated_'\$fnum'.card'

cp '$card/t2hkk_oscillated.card' '$card/t2hkk_oscillated_'\$fnum'.card' 

$Exec $card/t2hkk_oscillated_\$fnum'.card'
EOF

echo qsub -cwd $scriptname





#for j in {0..19} ; do
#	cp '$card/t2hkk_oscillated.card' $card/t2hkk_oscillated_$j'.card' 
#	first=$(expr $chunk \* $j)
#	last=$(expr $first + $chunk - 1)
#	echo $first $last 
#	sed -i 's:first.*:first $first:' $card/t2hkk_oscillated_$j'.card'
#	sed -i 's:last.*:last $last:' $card/t2hkk_oscillated_$j".card"
#
#	cat > $output/Roscillator_$j".sh" <<EOF
##!/bin/bash
#source /home/tboschi/OscAna/Setup_OscAna.sh
#$Oscillator $card/t2hkk_oscillated_$j.card
#EOF
#	echo -r Oscillator$j -o $output/Roscillator_$j.log -e $output/Roscillator_$j.err $output/Roscillator_$j".sh" >> /home/tboschi/jobManager/atmpd.list
#done
#
#if [ $rest -ne 0 ] ; then
#	cp "$card/t2hkk_oscillated.card" $card/t2hkk_oscillated_20.card 
#	first=$(expr $last + 1)
#	last=$(expr $first + $rest)
#	echo $first $last 
#	sed -i "s:first.*:first $first:" $card/t2hkk_oscillated_20.card
#	sed -i "s:last.*:last $last:" $card/t2hkk_oscillated_20.card
#fi
#
##smart submit to the queue system
#queues=(atmpd ALL all lowe calib)
#count=0
#list=$JM/${queues[$count]}.list
#echo 'Submitting to the queue system'
#while [ -s $list ] ; do
#	wcleft=$(wc -l $list)
#	jobsinqueue=$(qstat ${queues[$count]} | grep run | awk '{print $1;}')
#	echo "   jobs left to sumbit : " ${wcleft%%/*}
#	echo "   jobs in " ${queues[$count]} " : " $jobsinqueue
#	$Queue $(expr $MAX_JOBS - $jobsinqueue)
#        count=$(expr $count + 1)
#	count=$(expr $count % ${#queues[@]})
#	echo "count " $count
#	mv $list $JM/${queues[$count]}.list
#	list=$JM/${queues[$count]}.list
#done

#mv $cardS/t2hkk_oscillated*.card $output/
