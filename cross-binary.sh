#! /bin/bash

# this script detects the architecture of the machine
# running on it and launches the right binary
# if the binary does not exist, a generic binary is used

arch=$(gcc -march=native -Q --help=target | grep march | cut -f3)

exe=$1
shift 

echo Running on $(hostname) with $arch
if [ -s $PWD/bin/arch/"$exe"_$arch ] ; then
	$PWD/bin/arch/"$exe"_$arch "$@"
elif [ -s $PWD/bin/arch/"$exe" ] ; then
	$PWD/bin/arch/"$exe" "$@"
elif [ -s $PWD/bin/"$exe" ] ; then
	$PWD/bin/"$exe" "$@"
else
	echo There is no binary called \""$exe"\"
fi
