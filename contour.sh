#! /bin/bash

usage="Usage: $0 point_path
			
Create X2 contours

  parameters
    point_dir    output directory of complete fit, usually ends with ...point_xxx/
"

if [ "$#" -ne 1 ] ; then
	echo This script requires one argument. Check usage with
	echo $0 -h
	exit 1
elif [ "$1" == "-h" ] ; then
	echo "$usage"
	exit 0
elif [ ! -s "$1" ] ; then
	echo Argument \"$1\" is not valid. Check usage with
	echo $0 -h
	exit 1
fi

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

	$$PWD/bin/buildcontours $cont/all_$name.root $cont/$name.root $sens/oscillation.card

	rm $cont/all_$name.root
done

echo "DONE"
