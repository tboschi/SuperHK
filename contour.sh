#! /bin/bash

Contour=/data/tboschi/HKsens/OscAna/SuperHK/bin/buildcontours

point=$(cat $1)
point=(${point})

root=${1%.*}
name=${root##*/}
sens=${1%/*}
cont=$sens/../contours

for p in "${point[@]}"
do
	echo contouring point $p

	mkdir -p $cont/$name

	hadd -f $cont/$name/all.root $sens/$name_$p/*.root

	$Contour $cont/$name/all.root $cont/$name/contour.root $sens/oscillation.card

	rm $cont/$name/all.root
done

echo "DONE"
