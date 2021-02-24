#! /bin/bash

usage="usage: $0 <file> [ <file2> ]
			
Create CPV exclusion curves using the output of a CPV scan created
with trisens.sh. The input <file> is found in the output directory, 
with extension \".info\" and contains a list of the true deltaCP
points fitted.
Using the parameter names of trisens.sh, the file should look like

	<root>/<mh1>_<mh2>/sensitivity/CPV_scan.info

The output is a text file which will be saved in

	<root>/<mh1>_<mh2>/contours/CPV_scan.dat

If a second file <file2> is passed, then the script will create
the CPV curves for the other fit as well.
The second file should be a CPV scan made with a different fit mass
hierarchy, but same true mass hierarchy (see trisens.sh documentation).
The script will also create CPV exclusion curves for unknown mass
hierarchy of the common true mass hierarchy, and it will be saved in

	<root>/<mh1>/contours/CPV_scan.dat

"

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ] ; then
	echo This script requires one or two arguments. Check usage with
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

point=$(cat $1)
point=(${point})

root=${1%.*}
name=${root##*/}
cont=${root%/*}/../contours
root=$root"_"

point=("${point[@]/#/$root}")

mkdir -p $cont/
echo Creating exclusion for $cont

$PWD/bin/exclusion $cont/"$name".dat "${point[@]}"

# this is executed only if two info points are passed
if [ "$#" -eq 2 ] ; then
	$0 $2

	comm=${2%_*/sensitivity*}/contours
	cont2=${2%/*}/../contours

	mkdir -p $comm

	echo Creating unknown mass hierarchy exclusion for $comm
	$PWD/bin/unknown_mh $cont/"$name".dat $cont2/"$name".dat $comm/"$name".dat
fi
