#! /bin/bash

Chi2=/data/tboschi/HKsens/OscAna/SuperHK/bin/dropchi2
Sens=/data/tboschi/HKsens/OscAna/SuperHK/bin/dropsens
root=/data/tboschi/HKsenHK

#samples=(T2HK)

pp=false
while getopts 'pf:r:n:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		f) point="${OPTARG}" ;;
		n) name="${OPTARG}" ;;
		p) pp=true ;;
		*) exit 1 ;;
	esac
done

echo $root, $name
shift $((OPTIND-1))

root=${root%/}

if [ "$pp" = true ] ; then
	tgt="gaussian"
else
	tgt="uniform"
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

echo "creating dir" $name
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
