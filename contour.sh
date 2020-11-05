#! /bin/bash

usage="usage: $0 <point>
			
Create X2 contours for a fit result obtained with trisens.sh.
The input <point> is the path to the output folder. Using the
parameter names of trisens.sh, the output folder names should
look like

	<root>/<mh1>_<mh2>/sensitivity/point_<value>

and <value> is the number of the fitted point.
The X2 contours are saved as histograms in the ROOT file

	<root>/<mh1>_<mh2>/contours/point_<value>.root

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
