#! /bin/bash

Contour=/home/tboschi/OscAna/Osc3++/processing/build.contours/BuildContourPlots
root=/home/tboschi/data/
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

list=$(ls -d $root/sensitivity/*_*/)

for ff in "${samples[@]}" ; do
	echo Processing $ff set
	for ll in "${list[@]}" ; do
		dir=${ll%/}
		dir=${dir##*/}
		mkdir -p $root/contours/$dir

		#contour of pure files
		hadd -f $root/contours/$dir/all.$ff.root $root/sensitivity/$dir/SpaghettiSens.$ff.*.root
		$Contour $root/contours/$dir/all.$ff.root >& /dev/null
		mv ChiSquared.root $root/contours/$dir/uniform.$ff.root
		rm $root/contours/$dir/all.$ff.root

		#contour of penalised files
		hadd -f $root/contours/$dir/all.$ff.root $root/sensitivity/$dir/SpaghettiSens_penalised.$ff.*.root
		$Contour $root/contours/$dir/all.$ff.root >& /dev/null
		mv ChiSquared.root $root/contours/$dir/gaussian.$ff.root
		rm $root/contours/$dir/all.$ff.root
	done
done

echo "DONE"
