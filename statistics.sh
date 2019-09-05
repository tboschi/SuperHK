#! /bin/bash

#Osc=/home/tboschi/OscAna/hk.atm+beam/submitOsc.sh

Stat=/home/tboschi/OscAna/hk.atm+beam/submitSensStatOnly.sh
JM=/home/tboschi/jobManager/
Queue=$JM/submitToQueue.sh

root=/home/tboschi/data/
card=/home/tboschi/OscAna/hk.atm+beam/
MH=NH

while getopts 'inr:c:' flag; do
	case "${flag}" in
		i) MH=IH ;;
		n) MH=NH ;;
		r) root="${OPTARG}" ;;
		c) card="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

NFILES=$(ls $global/ | wc -l)


cp $card temp.card

#get ready to sensitivity

cp temp.card temp2.card
mv temp.card temp1.card

sed -i "s:InputToFit.*:InputToFit $root/global/NH/file.*.root:" temp2.card
if [[ $MH = "NH" ]]; then
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" temp2.card
else
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" temp2.card
fi

#perform fit
for i in {0..200}
do
	echo test point $i
	echo $Stat $NFILES temp1.card temp2.card $i $root/statistics/$MH/syst_t$i
done

#smart submit to the queue system
echo 'Submitting to the queue system'
queues=(atmpd ALL all lowe calib)
count=0
list=$JM/${queues[$count]}.list
echo $list
wc -l $list
#while $(wc -l < $list) ; do
while [ -s $list ] ; do
	#qstat -a | grep 
	echo $Queue
        count=$(expr $count + 1)
	echo mv $list $JM/${queues[$count]}.list
	list=$JM/${queues[$count]}.list
	echo $list $count
done

#wait until jobs are finished
while true ; do
	check=$(qstat -u tboschi | grep tboschi)
	if [[ $check ]]; then
		echo 'waiting 1min...'
		sleep 60
	else
		break
	fi
done

echo ./applypenalty cards/penalty.card

listSK=$(ls $root/statistics/$MH/syst_*/Spaghetti.SK.*.root) 
listT2HK=$(ls $root/statistics/$MH/syst_*/Spaghetti.T2HK.*.root) 
listSKT2HK=$(ls $root/statistics/$MH/syst_*/Spaghetti.SKT2HK.*.root) 

echo ./buildmulti.sh list_SK
echo ./buildmulti.sh list_T2HK
echo ./buildmulti.sh list_SKT2HK

echo ./exclusion list_SK
echo ./exclusion list_T2HK
echo ./exclusion list_SKT2HK

echo "DONE"
