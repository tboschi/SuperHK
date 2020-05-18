#! /bin/bash

Exclude=./bin/exclusion
root=/data/tboschi/HKsens/
samples=(T2HK)

aa=false
pp=false
tr=false
while getopts 'fr:pa' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		p) pp=true ;;
		a) aa=true ;;
		f) tr=true ;;
		*) exit 1 ;;
	esac
done

#root contains NH_NH
root=${root%/}

if [ "$tr" = true ] ; then
	Exclude=./bin/exclusion_filter
	tgt='gaussian_filter'
else
	Exclude=./bin/exclusion
	tgt='gaussian'
fi

mkdir -p $root/exclusion/

for ff in "${samples[@]}" ; do
	if [ "$pp" = true ] ; then
		$Exclude $root/sensitivity 'SpaghettiSens_penalised'.$ff
		mv Exclusion.dat $root/exclusion/$tgt.$ff.dat
	elif [ "$aa" = true ] ; then
		$Exclude $root/sensitivity 'SpaghettiSens_atmo'.$ff
		mv Exclusion.dat $root/exclusion/combined.$ff.dat
	else
		$Exclude $root/sensitivity 'SpaghettiSens'.$ff 
		mv Exclusion.dat $root/exclusion/uniform.$ff.dat
	fi
done

echo "DONE"
