#! /bin/bash

usage="usage: $0

Compile the fitter and the atmo_input binaries on some nodes on the cluster.
The script runs through all the nodes and compiles a new version of the binary
only if a new architecture is found. If it was not possible to ssh into any
node, then a generic version of the two binaries will be compiled, otherwise
a highly optimized version is creted.

Works with HTCondor or Slurm and requires shared file system between nodes.
"


if [ "$#" -gt 0 ] ; then
	echo No arguments required.
	echo "$usage" >&2
	exit 1
fi


# tool to determine arch
mkdir -p bin/arch
cat > bin/arch/arch.sh  << EOF
#! /bin/bash
hostname
gcc -march=native -Q --help=target | grep march
sleep 10
EOF

if condor_q &> /dev/null ; then
	hpc=$(condor_status | awk '/@/ {sub(/.*@/, "", $1); print $1}' | sort -u)
elif squeue &> /dev/null ; then
	hpc=$(sinfo -h -N -p nms_research,shared -o "%n" | sort -u)
else
	echo There is neither HTCondor nor Slurm on this machine. I am sorry, I cannot help you
	exit 1
fi

hpc=(${hpc})

echo Hosts found in cluster:
echo "    " "${hpc[@]}"
begin=$(date +%s)

for host in "${hpc[@]}"
do
	ssh -o UserKnownHostsFile=/dev/null \
	    -o LogLevel=ERROR \
	    -o PasswordAuthentication=no \
	    -o StrictHostKeyChecking=no \
	    $host /bin/bash << EOF
echo on \$(hostname)
if [ -f ~/.bash_profile ] ; then
	. ~/.bash_profile
elif [ -f ~/.profile ] ; then
	. ~/.profile
fi
cd $PWD
arch=\$(gcc -march=native -Q --help=target | grep march | cut -f3)
tgt=bin/arch/fitter_\$arch
edit=$(date +%s)
if [ -s \$tgt ] && [ \$(date +%s -r \$tgt) -gt $begin ] ; then
	exit
fi
echo compiling for \$arch
make clean
make APP=fitter
mv bin/fitter \$tgt
make APP=atmo_input
mv bin/atmo_input bin/arch/atmo_input_\$arch
EOF
done

echo on localhost
make clean

# compile fitter and atmo_input with generic architecture
# if it wasn't possible to compile on the cluster
if ls bin/arch/fitter_* &> /dev/null ; then
	make APP=fitter ARCH=
fi
if ls bin/arch/atmo_input_* &> /dev/null ; then
	make APP=atmo_input ARCH=
fi

#compile the rest
make
