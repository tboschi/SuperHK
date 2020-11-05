#! /bin/bash

# this script detects the architecture of the machine
# running on it and launches the right binary
# if the binary does not exist, a generic binary is used

arch=$(gcc -march=native -Q --help=target | grep march | cut -f3)

exe=$1
shift 

echo Running on $arch
if [ -s $PWD/bin/arch/"$exe"_$arch ] ; then
	$PWD/bin/arch/"$exe"_$arch "$@"
else
	$PWD/bin/"$exe" "$@"
fi
