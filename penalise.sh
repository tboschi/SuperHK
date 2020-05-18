#! /bin/bash

Penalty=/data/tboschi/HKsens/OscAna/SuperHK/bin/addpenalty
cardS=/data/tboschi/HKsens/OscAna/SuperHK/cards
samples=(T2HK)
csub=condor_submit

while getopts 'r:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

shift $((OPTIND-1))

root=${root%/}
base=${root%/?H_?H}
base=${base##*/}
root=$root/sensitivity
cp $cardS/penalty_$base.card $root/penalty_x2.card

echo $root
echo $base

for point in "$@"
do
	list=("${list[@]}" "$root/$point")
done

if [ ${#list[@]} -eq 0 ]; then
	list=($(ls -d $root/*_*/))
fi

echo $list

for ff in "${samples[@]}" ; do
	echo Processing $ff set
	for ll in "${list[@]}" ; do
		dir=${ll%/}
		dir=${dir##*/}

		nfiles=$(find $root/$dir/ -name "SpaghettiSens.$ff.*.root" | wc -l)
		if [ $nfiles -eq '0' ] ; then
			continue
		fi

		tupdt=1
		tlast=0
		npenals=$(find $root/$dir/ -name "SpaghettiSens_penalised.$ff.*.root" | wc -l)
		if [ $npenals -gt 0 ] ; then
			last=$(ls -t $root/$dir/SpaghettiSens_penalised.$ff.*.root | head -n1)
			updt=$(ls -t $root/$dir/SpaghettiSens.$ff.*.root | head -n1)
			tlast=$(date -r $last +%s)
			tupdt=$(date -r $updt +%s)
		fi
		echo $tlast , $tupdt

		if [ $nfiles -ne $npenals -o $tupdt -gt $tlast ] ; then
			echo 'Penalising' $dir

			#penalise
			rm -f $root/$dir/'SpaghettiSens_penalised.'$ff'.*.root'
			echo $scriptname
			scriptname=$root/$dir/Rpenalise_$ff'.sub'
			cat > $scriptname << EOF
# script submission for condor
# sumbit with --
#	condor_submit $scriptname


executable		= $Penalty
arguments		= $root/penalty_x2.card $root/$dir/
getenv			= True
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $root/$dir/Lpenalise_$ff.log
error			= $root/$dir/Lpenalise_$ff.log

queue

EOF

#'SpaghettiSens.$ff.*.root
			$csub $scriptname
			#$Penalty $root/penalty_x2.card $root/$dir/SpaghettiSens.$ff.*.root
		else
			echo 'Skipping' $dir
		fi
	done
done

echo 'DONE'
