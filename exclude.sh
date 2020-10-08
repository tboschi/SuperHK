#! /bin/bash

Exclude=/data/tboschi/HKsens/OscAna/SuperHK/bin/exclusion

point=$(cat $1)
point=(${point})


root=${1%.*}
name=${root##*/}
cont=${root%/*}/../contours
root=$root"_"

point=("${point[@]/#/$root}")

mkdir -p $cont/

$Exclude $cont/$name.dat "${point[@]}"

echo DONE
