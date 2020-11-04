#! /bin/bash

usage="Usage: $0 info_file [ info_file2 ]
			
Create CPV exclusion curves

  parameters
    info_file    path to file with .info extension in output folder
    info_file2   path to file with .info extension and different fit mass hierarchy,
    		 but with the same true mass hierarchy; exclusion lines with unknown
		 mass hierarchy will be created as well under a new folder
"

if [ "$#" -lt 1 ] ; then
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
