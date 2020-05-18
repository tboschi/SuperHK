#! /bin/bash

get=./bin/getepsilon
root=/data/tboschi/HKsens/

while getopts 'r:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

#root contains NH_NH
root=${root%/}
out=$root/contours/point_49773/errors.dat
root=$root/sensitivity/point_49773/SpaghettiSens_penalised.T2HK.\*.root


echo $get "$root" "$out"
$get "$root" "$out"

echo "DONE"
