#! /bin/bash

Contour=/data/tboschi/HKsens/OscAna/SuperHK/bin/buildcontours

full=${1%/}
name=${full##*/}
point=${name##*_}
sens=${full%/*}
cont=$sens/../contours

for p in "${point[@]}"
do
	echo contouring point $p

	mkdir -p $cont/

	hadd -f $cont/all_$name.root $sens/$name/*.root

	$Contour $cont/all_$name.root $cont/$name.root $sens/oscillation.card

	rm $cont/all_$name.root
done

echo "DONE"
