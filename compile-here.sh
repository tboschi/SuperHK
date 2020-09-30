#! /bin/bash


for host in "$@"
do
	ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no $host /bin/bash << EOF
echo on \$(hostname)
source .profile
cd $PWD
arch=\$(gcc -march=native -Q --help=target | grep march | cut -f3)
echo compiling for \$arch
make clean
make APP=fitter
mv bin/fitter bin/arch/fitter_\$arch
EOF
done

echo on localhost
make clean
make
