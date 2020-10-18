#! /bin/bash

arch=$(gcc -march=native -Q --help=target | grep march | cut -f3)

echo Running on $arch
$PWD/bin/arch/atmo_input_$arch "$@"
