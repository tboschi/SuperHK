#! /bin/bash

# this script detects the architecture of the machine
# running on it and launches the right binary
# if the binary does not exist, a generic binary is used

arch=$(gcc -march=native -Q --help=target | grep march | cut -f3)

echo Running on $arch
if [ -s $PWD/bin/arch/atmo_input_$arch ] ; then
	$PWD/bin/arch/atmo_input_$arch "$@"
else
	$PWD/bin/atmo_input "$@"
