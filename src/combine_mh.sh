#! /bin/bash

Unknown=$PWD/bin/unknown_mh

root=${1%/}

nh="$root"_NH/contours/CPV_scan.dat
ih="$root"_IH/contours/CPV_scan.dat

if [ ! -s $nh ]; then 
	echo No file $nh
	exit 1
fi

if [ ! -s $ih ]; then 
	echo No file $ih
	exit 1
fi

mkdir -p $root
out="$root"/CPV_unknown.dat

$Unknown $nh $ih $out
