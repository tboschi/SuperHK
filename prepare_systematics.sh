#! /bin/bash

while getopts 'r:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

shift $((OPTIND-1))

root=${root%/}
root=$root/systematics
echo $root

fhc1Re_i=$root/1Re.*.t2k.root
fhc1Rm_i=$root/1Rmu.*.t2k.root
rhc1Re_i=$root/RHC1Re.*.t2k.root
rhc1Rm_i=$root/RHC1Rmu.*.t2k.root

fhc1Re_o=$root/FHC1Re.fij.t2k_spline.root
fhc1Rm_o=$root/FHC1Rmu.fij.t2k_spline.root
rhc1Re_o=$root/RHC1Re.fij.t2k_spline.root
rhc1Rm_o=$root/RHC1Rmu.fij.t2k_spline.root

#banff=("banff" "Postfit")
#for nn in "${banff[@]}" ; do
#	matrixBanff=$root/systematics/*$nn*.root
#	echo banff $nn $matrixBanff
#	if [ -e $matrixBanff ] ; then
#		echo $matrixBanff found
#		break
#	fi
#done

#skdet=("SKJoint" "skdet")
#for nn in "${skdet[@]}" ; do
#	matrixSKdet=$root/systematics/*$nn*.root
#	echo skdet $nn $matrixSKdet
#	if [ -e $matrixSKdet ] ; then
#		echo $matrixSKdet found
#		break
#	fi
#done

./bin/purifysystematic $fhc1Re_i $fhc1Re_o
./bin/purifysystematic $fhc1Rm_i $fhc1Rm_o
./bin/purifysystematic $rhc1Re_i $rhc1Re_o
./bin/purifysystematic $rhc1Rm_i $rhc1Rm_o


matrix=()
for f in "$@"
do
	echo $root/$f
	matrix=("${matrix[@]}" "$root/$f")
	echo ${matrix[@]}
done

if [ ${#matrix[@]} -eq 0 ]; then
	echo "Nothing else to do"
else
	./bin/addmatrix ${matrix[@]}
	mv combinedmatrix.root $root
fi
