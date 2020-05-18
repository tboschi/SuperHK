#! /bin/bash

#run GlobalOsc efficiently

Exec=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalOsc

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
global=$global/$MH/global
card=$global/global.card

#name=$(grep OutputOscFile $card | cut -d'"' -f2)
#output=${name%.*}

#echo sed -i "s:OutputOscFile.*:OutputOscFile \"$output.root\":" $card

if [[ $MH = "NH" ]]; then
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" $card
else
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" $card
fi

nameExec=${Exec##*Global}
rm -f $global/L$nameExec*log
scriptname=$global/'R'$nameExec'.sh'
MAX_JOBS=300

cat > $scriptname << EOF
#$ -S /bin/bash
#$ -N L$nameExec
#$ -t 1-$NFILES
#$ -tc $MAX_JOBS
#$ -o $global/\$JOB_NAME.\$TASK_ID.log
#$ -j y

$Exec \$(expr \$SGE_TASK_ID - 1) \$(expr \$SGE_TASK_LAST - \$SGE_TASK_FIRST) $global $card
EOF

echo qsub -cwd $scriptname

#$Osc $card $global/global $NFILES

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

echo "global finished"
