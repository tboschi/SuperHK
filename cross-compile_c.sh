#! /bin/bash


# tool to determine arch
mkdir -p bin/arch
cat > bin/arch/arch.sh  << EOF
#! /bin/bash
hostname
gcc -march=native -Q --help=target | grep march
sleep 10
EOF

info=condor_status
hpc=$($info | awk -v m=$tag '/ph\.qmul\.ac\.uk/ {sub(/.*@/, ""); print $1} ' | sort -u)
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
EOF
done

echo on localhost
make clean
make
