#! /bin/bash


mkdir -p bin/arch
cat > bin/arch/arch.sh  << EOF
#! /bin/bash
hostname
gcc -march=native -Q --help=target | grep march
sleep 10
EOF

cat > bin/arch/arch.sub << EOF
#! /bin/bash
#SBATCH --array=0-111
#SBATCH --job-name=print_arc
#SBATCH -o bin/arch/arch_list.%a
#SBATCH -p nms_research,shared
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1

srun bin/arch/arch.sh
EOF

chmod u+x bin/arch/arch.sh bin/arch/arch.sub

sbatch bin/arch/arch.sub

while [ $(squeue -u $USER | wc -l) -gt 1 ] ; do
	sleep 1
done

arches=$(grep -h march bin/arch/arch_list.* | cut -f3 | sort -u)
arches=(${arches})

echo Architectures detected in the cluster: $arches

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
