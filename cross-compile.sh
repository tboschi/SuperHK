#! /bin/bash

record="bin/arch/nodes_arches"
usage="usage: $0 [<options>]

Compile the fitter and the atmo_input binaries on some nodes on the cluster.
The script runs through all the nodes and compiles a new version of the binary
only if a new architecture is found. If it was not possible to ssh into any
node, then a generic version of the two binaries will be compiled, otherwise
a highly optimized version is creted.
After running the first time, the script will create the minimal list of nodes
for a thorough compilation in
	$PWD/$record
and ssh only into those nodes in subsequent compilations.

Works with HTCondor or Slurm and requires shared file system between nodes.

Optional parameters
   -a			 visit all nodes on the cluster. This updates the list
   			 of nodes and architectures as well.

   -x <node>[,<nodes>]	 specify a node or list of nodes to exclude
   			 from compilation. For a list use commas and
			 not white spaces. Wildcards are supported.

   -h			 print this message and quit.

"

all=false
exclude=""
while getopts 'ax:h' flag; do
	case "${flag}" in
		a) all=true ;;
		x) exclude="${OPTARG}" ;;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done

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

if [ "$all" == "true" ] ; then
	rm -r $record
fi

if [ -s $record ] ; then
	hpc=$(cat $record | cut -f1 -d' ')
else
	touch $record
fi

hpc=(${hpc})
if [ -n "$exclude" ]; then
	exclude=(${exclude//,/ })
	for e in "${exclude[@]}" ; do
		hpc=("${hpc[@]/$e}")
	done
fi

echo Beginning compilation on the following nodes:
echo "    " "${hpc[@]}"
begin=$(date +%s)

if [ -n "$visit" ]; then
	hpc=(${visit//,/ })
	echo "Visiting:"
	echo "    " "${hpc[@]}"
fi


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
if ! grep -q "\$arch" $record ; then
	echo \$(hostname)  "\$arch" >> $record
fi
echo compiling for \$arch
make clean
make -j4 APP=fitter
mv bin/fitter \$tgt
make -j4 APP=atmo_input
mv bin/atmo_input bin/arch/atmo_input_\$arch
EOF
done

echo on localhost

# compile fitter and atmo_input with generic architecture
# in case it wasn't possible to compile them on the cluster
make clean
make -j4 APP=fitter ARCH=
mv bin/fitter bin/arch/fitter

make clean
make -j4 APP=atmo_input ARCH=
mv bin/fitter bin/arch/fitter

#compile the rest
make clean
make -j4
