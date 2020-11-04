#! /bin/bash

arch=$(gcc -march=native -Q --help=target | grep march | cut -f3)

echo Running on $arch
if [ -s $PWD/bin/arch/fitter_$arch ] ; then
	$PWD/bin/arch/fitter_$arch "$@"
else
	$PWD/bin/fitter	"$@"
fi
