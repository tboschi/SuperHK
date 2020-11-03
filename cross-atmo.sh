#! /bin/bash

arch=$(gcc -march=native -Q --help=target | grep march | cut -f3)

echo Running on $arch
if [ -s $PWD/bin/arch/atmo_input_$arch ] ; then
	$PWD/bin/arch/atmo_input_$arch "$@"
else
	$PWD/bin/atmo_input "$@"
