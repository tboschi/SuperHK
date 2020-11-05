#! /bin/bash

usage="usage: $0 -r root_folder -1 [NH | IH] [-N num_jobs]
              [-t fraction] [-v verbosity] [-h]

Pre compute the atmospheric sample, using HTCondor or Slurm.

  parameters
    -r root_folder  output folder
    -1 [NH | IH]    mass hierarchy of the observed/true data (NH or IH)

  optional parameters
    -N num_jobs     number of jobs to submit, default 360
    -t fraction	    specify fraction of data to fit as value [0,1], default 1 (full data)
    -v verbosity    specify a verbosity value
    -h		    show usage"

SCHED=""
if condor_q &> /dev/null ; then
	SCHED="HTCONDOR"
elif squeue &> /dev/null ; then
	SCHED="SLURM"
else
	echo There is neither HTCondor nor Slurm on this machine. I am sorry, I cannot help you
	exit 1
fi


Atmo=$PWD/cross-binary.sh
nameExec=atmo_input

card=$PWD/cards/multi.card
oscc=$PWD/cards/oscillation.card
atmo=$PWD/cards/atmo_sample.card

MAX_JOBS=500
MAX_QUEUE=500

root=""
data=""
#global=/data/tboschi
MH_1=""
verb="1"

while getopts 'r:1:N:t:v:h' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
		N) NJOBS="${OPTARG}" ;;
		t) stats="${OPTARG}" ;;
		v) verb="${OPTARG}" ;;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done

rm -f .reconstruction_files
rm -f .production_files

#first create output folders

if [ -z $root ] ; then
	echo You must specify a root folder with -r
	exit
fi

if [ -z $MH_1 ] ; then
	echo You must define mass hirerachy with -1 and -2
	exit
fi

# add PWD 
if [[ "$root" != /* ]] ; then
	root=$PWD/${root%/}
else
	root=${root%/}
fi


#define mass hierarchy to fit
root=$root/$MH_1
mkdir -p $root
rm -f $root/*.*

# copy cards to output folder
cp $card $oscc $atmo $root/
card=$root/${card##*/}
oscc=$root/${oscc##*/}
atmo=$root/${atmo##*/}

if [ -n $verb ] ; then
	echo Setting verbosity to $verb
	sed -i "s:verbose.*:verbose\t$verb:" $card $oscc $atmo
fi

if [ -z $NJOBS ] ; then
	NJOBS=$MAX_JOBS
fi

#use correct hierachies
if [ $MH_1 = "NH" ] ; then
	sed -i "s:^true_hierarchy.*:true_hierarchy\t\"normal\":" $card
elif [ $MH_1 = "IH" ] ; then
	sed -i "s:^true_hierarchy.*:true_hierarchy\t\"inverted\":" $card
fi

sed -i "/^beam_parameters/s:^:#:" $card
sed -i "/^fit_parameters/s:^:#:" $card

sed -i "/^#atmo_parameters/s:^#::" $card
sed -i "s:^atmo_parameters.*:atmo_parameters\t\"$atmo\":" $card
sed -i "s:^oscillation_parameters.*:oscillation_parameters\t\"$oscc\":" $card

sed -i "s:^output.*:output\t\"$root/atmo.root\":" $card



#update statistics
if [ -z $stats ] ; then
	sed -i "/^stats/s:^:#:" $atmo
else
	sed -i "/^#stats/s:^#::" $atmo
	sed -i "s:stats.*:stats\t0.$stats:" $atmo
fi


#atmo systematics

sys_atmo=$root/../../atmo_fij.root
sed -i "s:^systematic_file.*:systematic_file\t\"$sys_atmo\":" $atmo
sed -i "s:^systematic_tree.*:systematic_tree\t\"sigmatree\":" $atmo
#Atmo systematics
sed -i "/^#stats_only/s:^#::" $atmo

reco_atmo=$root'/../../*.sk4.*.root'
#MC inputs
sed -i "s:^MC_input.*:MC_input\t\"$reco_atmo\":"	$atmo
sed -i "s:^MC_tree_name.*:MC_tree_name\t\"osc_tuple\":"	$atmo

dens=$PWD'/data/PREM_25pts.dat'
prod=$PWD'/data/prod_honda/kam-ally-aa-*.d'
sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$atmo
sed -i "s:honda_production.*:honda_production\t\"$prod\":"	$atmo
sed -i "s:production_height.*:production_height\t15:"		$atmo


scriptname=$root/R$nameExec.sub

if [ "$SCHED" == "HTCONDOR" ] ; then
	sub=condor_submit
	cat > $scriptname << EOF
# script submission for condor
# sumbit with --
#	$sub $scriptname

executable		= $Atmo
arguments		= atmo_input \$(Process) $NJOBS $card
getenv			= True
#requirements		= HasAVXInstructions
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $root/L$nameExec.\$(Process).log
error			= $root/L$nameExec.\$(Process).log

queue $NJOBS

EOF
	echo Launching $NJOBS jobs with HTCondor
elif [ "$SCHED" == "SLURM" ] ; then
	sub=sbatch
	cat > $scriptname << EOF
#! /bin/bash
# script submission for SLURM
# sumbit with --
#	$sub $scriptname

#SBATCH --array=0-$((NJOBS - 1))
#SBATCH --job-name=$nameExec
#SBATCH -o $root/L$nameExec.%a.log
#SBATCH -p nms_research,shared
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1

srun $Atmo atmo_input \$SLURM_ARRAY_TASK_ID $NJOBS $card

EOF
	echo Launching $NJOBS jobs with Slurm
fi

echo Submitting jobs
$sub $scriptname
