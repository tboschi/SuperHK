#! /bin/bash


mkdir -p bin/arch
cat > bin/arch/arch.sh  << EOF
#! /bin/bash
hostname
gcc -march=native -Q --help=target | grep march
sleep 10
EOF

cat > bin/arch/arch.sub << EOF
executable		= bin/arch/arch.sh
output			= bin/arch/arch_list.\$(Process)
getenv			= True
initialdir		= $PWD
queue 1000
EOF

chmod u+x bin/arch/arch.sh bin/arch/arch.sub

condor_submit bin/arch/arch.sub

while [ $(condor_q -all -format "%d\n" ProcId |  wc -l) -gt 0 ]; do
	sleep 1
	echo waiting...
done

arches=$(grep -h march bin/arch/arch_list.* | cut -f3 | sort -u)
arches=(${arches})

echo Architectures detected in the cluster: ${arches[@]}

for arch in "${arches[@]}" ; do

	host=$(grep -h $arch -B1 bin/arch/arch_list.* | head -n1)
	echo compiling for $arch on $host

	ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no $host /bin/bash << EOF
source .bash_profile
cd $PWD
make clean
make APP=fitter
mv bin/fitter bin/arch/fitter_$arch
EOF
done
