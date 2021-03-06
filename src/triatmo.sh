#! /bin/bash

usage="
usage: $0 -r <root> -1 <mh> [<options>]
              [-t fraction] [-v verbosity] [-h]

Pre-compute the atmospheric sample, using HTCondor or Slurm.
The <root> folder is the location of \"reconstruction_atmo\"
in the main analysis folder and <mh> is the desired mass
hierarchy.

Optional parameters
    -N <jobs>       number of jobs to submit to the cluster; the default value
                    is 360. There will be <jobs> output files in the end.
    -t <stat>	    specify fraction of data to fit as float value between
    -w <dir>        specify a different directory for log files as some file shared
    		    systems redirect output differently
    -v <verb>       specify a verbosity value where <verb> is an integer number;
    		    the greater the number, the higher the verbosity.
    -h		    print this message.
"

SCHED=""
if condor_q &> /dev/null ; then
	SCHED="HTCONDOR"
elif squeue &> /dev/null ; then
	SCHED="SLURM"
elif qstat -Q &> /dev/null ; then
    SCHED="PBS"
elif qstat -sql &> /dev/null ; then
    SCHED="SGE"
else
	echo "There is neither HTCondor nor Slurm nor qsub \(SGE or PBS\) on this machine. \
		Iam sorry, I cannot help you"
	exit 1
fi


Bin=$PWD/bin/cross-binary
Atmo=atmo_input

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
logr=""

while getopts 'r:1:w:N:t:v:h' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
		w) logr="${OPTARG}" ;;
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
root=$root/pre/$MH_1
mkdir -p $root
rm -f $root/*.*

if [ -n "$logr" ] ; then
	logr=$logr/pre/$MH_1/sensitivity
	mkdir -p $logr
	echo Log files are stored in $logr
else
	logr=$root
fi


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
# block pre computed entries
sed -i "/^pre_/s:^:#:" $atmo

dens=$PWD'/data/PREM_25pts.dat'
prod=$PWD'/data/prod_honda/kam-ally-aa-*.d'
sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$atmo
sed -i "s:honda_production.*:honda_production\t\"$prod\":"	$atmo
sed -i "s:production_height.*:production_height\t15:"		$atmo


scriptname=$root/R$Atmo.sub

case "$SCHED" in 
	"HTCONDOR")
		sub=condor_submit
		generate=$PWD/src/htcondor.submit
		;;
	"SLURM")
		sub=sbatch
		generate=$PWD/src/slurm.submit
		;;
	"PBS")
		sub=qsub
		generate=$PWD/src/pbs.submit
		;;
	"SGE")
		sub=qsub
		generate=$PWD/src/sge.submit
		;;
esac

$generate $Bin $Atmo $NJOBS $card $logr > $scriptname

echo Submitting $NJOBS jobs with $SCHED
$sub $scriptname


#if [ "$SCHED" == "HTCONDOR" ] ; then
#	sub=condor_submit
#	cat > $scriptname << EOF
## script submission for condor
## sumbit with --
##	$sub $scriptname
#
#executable		= $Atmo
#arguments		= atmo_input \$(Process) $NJOBS $card
#getenv			= True
##requirements		= HasAVXInstructions
#should_transfer_files	= IF_NEEDED
#when_to_transfer_output	= ON_EXIT
#initialdir		= $PWD
#output			= $logr/L$Atmo.\$(Process).log
#error			= $logr/L$Atmo.\$(Process).log
#
#queue $NJOBS
#
#EOF
#	echo Launching $NJOBS jobs with HTCondor
#elif [ "$SCHED" == "SLURM" ] ; then
#	sub=sbatch
#	cat > $scriptname << EOF
##! /bin/bash
## script submission for SLURM
## sumbit with --
##	$sub $scriptname
#
##SBATCH --array=0-$((NJOBS - 1))
##SBATCH --job-name=$Atmo
##SBATCH -o $logr/L$Atmo.%a.log
##SBATCH -p nms_research,shared
##SBATCH --time=3-0
##SBATCH --cpus-per-task=1
#
#srun $Atmo atmo_input \$SLURM_ARRAY_TASK_ID $NJOBS $card
#
#EOF
#	echo Launching $NJOBS jobs with Slurm
#elif [ "$SCHED" == "qsub" ] ; then
#    sub="qsub -t 1-$NJOBS -l sps=1 -o ${logr}/Logfile$Atmo.\$TASK_ID.log -e ${logr}/Logfile$Atmo.\$TASK_ID.log"
#    cat > $scriptname << EOF
#
##script submission for qsub
##submit with qsub -t 0-$NJOBS $scriptname
#cd /sps/t2k/lmuntean/HyperK/SuperHK/
#$Atmo atmo_input \$((SGE_TASK_ID-1)) $NJOBS $card 2>&1 | tee ${logr}/L$Atmo.\$((SGE_TASK_ID-1)).log
#
#EOF
#
#fi
