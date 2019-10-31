#! /bin/bash

Exclude=./bin/exclusion
root=/home/tboschi/data/
samples=(T2HK)

pp=false
while getopts 'r:p' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		p) pp=true ;;
		*) exit 1 ;;
	esac
done

#root contains NH_NH
root=${root%/}

mkdir -p $root/exclusion/

for ff in "${samples[@]}" ; do
	if [ "$pp" = true ] ; then
		$Exclude $root/sensitivity penalised
		mv Exclusion.dat $root/exclusion/gaussian.$ff.dat
	else
		$Exclude $root/sensitivity
		mv Exclusion.dat $root/exclusion/uniform.$ff.dat
	fi
done

echo "DONE"
