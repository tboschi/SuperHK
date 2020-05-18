#! /bin/bash

#oscillate files from GlobalOsc

Exec=/data/tboschi/HKsens/OscAna/SuperHK/bin/oscillator
card=/data/tboschi/HKsens/OscAna/SuperHK/cards
csub=condor_submit

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
reco=${global%/}/../../reconstruction
NFILES=$(ls $input/*.root | wc -l)
#NFILES=247

echo $input, $output, $reco
mkdir -p $output

cp $card/t2hkk_.card $output/oscillation.card
card=$output/oscillation.card

sed -i "s:input.*:input\t\"$input/*.root\":"	$card
sed -i "s:output.*:output\t\"$output\":"	$card
sed -i "s:pathtoreconstruction:$reco:"		$card

ih=$(grep InvertedHierarchy $input/global.card | cut -d ' ' -f 2)

if [ $ih -eq '0' ]; then
	sed -i "s:hierarchy.*:hierarchy\t\"normal\":"	$card
elif [ $ih -eq '1' ]; then
	sed -i "s:hierarchy.*:hierarchy\t\"inverted\":"	$card
fi

nameExec=${Exec##*bin/}
rm -f $output/R$nameExec*.sh

#oscillate beam part, updating cards with new systematics
N_JOBS=20
chunk=$(expr $NFILES / $N_JOBS)
resto=$(expr $NFILES % $N_JOBS)

echo $chunk, $resto

first=0
for id in $(seq 1 $N_JOBS) ; do
	if [ $id -le $resto ] ; then
		leg=$chunk
	else
		leg=$(expr $chunk - 1)
	fi

	last=$(expr $first + $leg)
	echo ID $id, $first --\> $last
	scriptname=$output/R$nameExec.$id.sub
	cat > $scriptname << EOF
# script submission for condor
# sumbit with --
#	condor_submit $scriptname


executable		= $Exec
arguments		= $first $last $card
getenv			= True
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $output/L$nameExec.$id.log
error			= $output/L$nameExec.$id.log

queue
EOF
	first=$(expr $last + 1)
	$csub $scriptname
done
echo "global finished"
