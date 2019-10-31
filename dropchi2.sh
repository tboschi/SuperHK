#! /bin/bash

Drop=/home/tboschi/OscAna/SuperHK/bin/dropchi2
samples=(T2HK)
root=/home/tboschi/data/

while getopts 'r:p:n:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		p) point="${OPTARG}" ;;
		n) name="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

shift $((OPTIND-1))

root=${root%/}

for ss in "$@"
do
	sets=("${sets[@]}" "$root/$ss/contours/$point/gaussian.T2HK.root")
done

echo "${sets[@]}"
$Drop ${sets[@]}

if [ $name ]
then
	mv X2minCP_all.dat $name'_X2minCP_all.dat'
	mv X2minM23_all.dat $name'_X2minM23_all.dat'
	mv X2minS13_all.dat $name'_X2minS13_all.dat'
	mv X2minS23_all.dat $name'_X2minS23_all.dat'
	mv X2minCP_diff.dat $name'_X2minCP_diff.dat'
	mv X2minM23_diff.dat $name'_X2minM23_diff.dat'
	mv X2minS13_diff.dat $name'_X2minS13_diff.dat'
	mv X2minS23_diff.dat $name'_X2minS23_diff.dat'
fi

echo "DONE"
