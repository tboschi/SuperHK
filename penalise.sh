#! /bin/bash

Penalty=/home/tboschi/OscAna/SuperHK/bin/addpenalty
cardS=/home/tboschi/OscAna/SuperHK/cards
samples=(T2HK)

while getopts 'r:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

root=${root%/}
base=${root%/?H_?H}
base=${base##*/}
cp $cardS/penalty_$base.card $root/sensitivity/penalty_sensitivity.card

list=$(ls -d $root/sensitivity/*_*/)

for ff in "${samples[@]}" ; do
	echo Processing $ff set
	for ll in "${list[@]}" ; do
		dir=${ll%/}
		dir=${dir##*/}
		#penalise
		rm -f $root/sensitivity/$dir/SpaghettiSens.$ff.*penalised*.root
		$Penalty $root/sensitivity/penalty_sensitivity.card $root/sensitivity/$dir/SpaghettiSens.$ff.*.root

	done
done

echo "DONE"
