#! /bin/bash

usage="Usage: $0

Compile the fitter binary on some nodes on the cluster.
The script runs through all the nodes and compiles a new version of the binary
only if a new architecture is found.

Works with HTCondor or Slurm and requires shared file system between nodes."


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
	hpc=$(condor_status | awk -v m=$tag '/ph\.qmul\.ac\.uk/ {sub(/.*@/, ""); print $1} ' | sort -u)
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
	    -o UserKnownHostsFile=/dev/null \
	    -o StrictHostKeyChecking=no \
	    $host /bin/bash << EOF
echo on \$(hostname)
source .profile
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
make
