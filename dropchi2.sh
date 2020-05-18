#! /bin/bash

Chi2=/data/tboschi/HKsens/OscAna/SuperHK/bin/dropchi2
Sens=/data/tboschi/HKsens/OscAna/SuperHK/bin/dropsens
root=/data/tboschi/HKsenHK

#samples=(T2HK)

ff=false
while getopts 'fr:p:n:' flag; do
	case "${flag}" in
		f) ff=true ;;
		r) root="${OPTARG}" ;;
		p) point="${OPTARG}" ;;
		n) name="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

shift $((OPTIND-1))

root=${root%/}

if [ "$ff" = true ]; then
	tgt='gaussian_filter'
else
	tgt='gaussian'
fi

for ss in "$@"
do
	setchi2=("${setchi2[@]}" "$root/$ss/contours/$point/$tgt.root")
	setsens=("${setsens[@]}" "$root/$ss/exclusion/$tgt.T2HK.dat")
done

if [ -s ${setchi2[0]} ] ; then
	echo "chi2"
	echo ${setchi2[@]}
	echo $Chi2 ${setchi2[@]}
	$Chi2 ${setchi2[@]}
fi

echo ""

if [ -s ${setsens[0]} ] ; then
	echo "sens"
	echo ${setsens[@]}
	echo $Sens ${setsens[@]}
	$Sens ${setsens[@]}
fi

mkdir -p plot/$name

if [ $name ]
then
	x2file=($(ls X2*.dat))

	for f in ${x2file[@]}; do
		echo $f
		mv $f plot/$name/$name'_'$f
	done

	ssfile=($(ls exclusion_*.dat))

	for f in ${ssfile[@]}; do
		echo $f
		mv $f plot/$name/$name'_'$f
	done
fi

echo "DONE"
