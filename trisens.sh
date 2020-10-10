#! /bin/bash

usage="Usage: $0 -r root_folder -1 [NH | IH] -2 [NH | IH] [-N num_jobs] [-s] 
              [-t fraction] [-f scan_mode] [-m matrix] [-v verbosity] [-h]

Submit fitter jobs to a scheduler, HTCondor or Slurm.

  parameters
    -r root_folder  main working folder, must have a \"systematics\" sub folder
    -1 [NH | IH]    mass hierarchy of the observed/true data (NH or IH)
    -2 [NH | IH]    mass hierarchy of the expected/to fit data (NH or IH)

  optional parameters
    -N num_jobs     number of jobs to submit, default 360
    -s		    does not fit systematic parameters, only statistics
    -t fraction	    specify fraction of data to fit as value [0,1], default 1 (full data)
    -f scan_mode    special scan mode, (CPV or octant)
    -m matrix	    specify matrix type for beam sample, default correlation
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



Sens=$PWD/cross-fitter.sh
nameExec=${Sens##*/}
nameExec=${nameExec%.*}

card=$PWD/cards/multi.card
fitc=$PWD/cards/fit_options.card
oscc=$PWD/cards/oscillation.card
beam=$PWD/cards/beam_sample.card
atmo=$PWD/cards/atmo_sample.card

MAX_JOBS=500
MAX_QUEUE=500

root=""
data=""
#global=/data/tboschi
MH_1=""
MH_2=""
ss=false	#false if systematic fit, true if only stat fit
ff=""		#empty if systematic fit, CPV if fast fit for dCP scan
mtype="correlation"
verb="1"

while getopts 'r:d:1:2:N:t:m:sf:v:h' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		2) MH_2="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
		d) data="${OPTARG}" ;;
		m) mtype="${OPTARG}" ;;
		N) NJOBS="${OPTARG}" ;;
		t) stats="${OPTARG}" ;;
		s) ss=true ;;
		f) ff="${OPTARG}" ;;
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

if [ -z $MH_1 ] || [ -z $MH_2 ] ; then
	echo You must define mass hirerachy for true and fitted samples with -1 and -2
	exit
fi

if [ -z $data ] ; then
	echo You must define a sample to fit with -d
	exit
fi

# add PWD 
if [[ "$root" != /* ]] ; then
	root=$PWD/${root%/}
else
	root=${root%/}
fi


#define mass hierarchy to fit
mhfit=$MH_1"_"$MH_2
root=$root/$mhfit

mkdir -p $root/sensitivity/

# copy cards to output folder
cp $card $fitc $oscc $beam $atmo $root/sensitivity/
card=$root/sensitivity/${card##*/}
fitc=$root/sensitivity/${fitc##*/}
oscc=$root/sensitivity/${oscc##*/}
beam=$root/sensitivity/${beam##*/}
atmo=$root/sensitivity/${atmo##*/}

if [ -n $verb ] ; then
	echo Setting verbosity to $verb
	sed -i "s:verbose.*:verbose\t$verb:" $card $fitc $oscc $beam $atmo
fi

case "$data" in
	"beam")
		echo beam
		sed -i "/^#beam_parameters/s:^#::" $card
		sed -i "/^atmo_parameters/s:^:#:" $card
		;;
	"atmo")
		echo atmo
		sed -i "/^#atmo_parameters/s:^#::" $card
		sed -i "/^beam_parameters/s:^:#:" $card
		;;
	"comb")
		echo comb
		sed -i "/^#beam_parameters/s:^#::" $card
		sed -i "/^#atmo_parameters/s:^#::" $card
		;;
esac

##defines input list
#input=$root/sensitivity/fit.list
#sed -i "s:^input.*:input\t\"$input\":" $card
#
### get number of files, atm there is one job per file with this definition
#if ls $global/$MH_1/global/*.root > $input ; then
#	NFILES=$(cat $input | wc -l)
#else
#	NFILES=300
#fi

if [ -z $NJOBS ] ; then
	NJOBS=360
fi


# decide if normal fit or fast fit for sensitivity
if [ -n "$ff" ] ; then
	pinfo=$ff"_scan.info"
	tname=$ff"_scan_"
	parm="CP"
	sed -i "/^#scan/s:^#::" $card
	sed -i "s:^scan.*:scan\t$ff:" $card
	sed -i "/^#point/s:^#::" $card
else
	pinfo="point.info"
	tname="point_"
	parm=""
	sed -i "/^scan/s:^:#:" $card
fi


#use correct hierachies
if [ $MH_1 = "NH" ] ; then
	sed -i "s:^true_hierarchy.*:true_hierarchy\t\"normal\":" $card
elif [ $MH_1 = "IH" ] ; then
	sed -i "s:^true_hierarchy.*:true_hierarchy\t\"inverted\":" $card
fi

if [ $MH_2 = "NH" ] ; then
	sed -i "s:^fit_hierarchy.*:fit_hierarchy\t\"normal\":" $card
elif [ $MH_2 = "IH" ] ; then
	sed -i "s:^fit_hierarchy.*:fit_hierarchy\t\"inverted\":" $card
fi



#update statistics
if [ -z $stats ] ; then
	sed -i "/^stats/s:^:#:" $beam $atmo
else
	sed -i "/^#stats/s:^#::" $beam $atmo
	sed -i "s:stats.*:stats\t0.$stats:" $beam $atmo
fi


#beam systematics

corr_beam=$root/../systematics/combinedmatrix.root
corr_atmo=$root/../systematics/atmo_corr.root
mtype="correlation"

# need to know correlation matrix to know number of systematics
sed -i "s:^corr_file.*:corr_file\t\"$corr_beam\":" $beam
sed -i "s:^corr_name.*:corr_name\t\"$mtype\":"  $beam

sys_E_FHC=$root/../systematics/FHC1Re.fij.t2k_spline.root
sys_E_RHC=$root/../systematics/RHC1Re.fij.t2k_spline.root
sys_M_FHC=$root/../systematics/FHC1Rmu.fij.t2k_spline.root
sys_M_RHC=$root/../systematics/RHC1Rmu.fij.t2k_spline.root
#just stats, comment systematics
if [ "$ss" = true ] ; then
	echo "statistics only fit"
	sed -i "/^systematic_/s:^:#:" $beam
else # update systematics
	sed -i "/^#systematic_/s:^#::" $beam
	sed -i "s:^systematic_E_FHC.*:systematic_E_FHC \"$sys_E_FHC\":" $beam
	sed -i "s:^systematic_E_RHC.*:systematic_E_RHC \"$sys_E_RHC\":" $beam
	sed -i "s:^systematic_M_FHC.*:systematic_M_FHC \"$sys_M_FHC\":" $beam
	sed -i "s:^systematic_M_RHC.*:systematic_M_RHC \"$sys_M_RHC\":" $beam
fi

reco_beam=$root'/../../reconstruction_beam/syst_*.card'
sed -i "s:^reco_input.*:reco_input\t\"$reco_beam\":" $beam

dens=$PWD'/data/DensityProfileTochibora.dat'
sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$beam


#atmo systematics

sys_atmo=$root/../systematics/atmo_fij.root
sed -i "s:^systematic_file.*:systematic_file\t\"$sys_atmo\":" $atmo
sed -i "s:^systematic_tree.*:systematic_tree\t\"sigmatree\":" $atmo
#Atmo systematics
if [ "$ss" = true ] ; then
	echo "statistics only fit"
	sed -i "/^#stats_only/s:^#::" $atmo
else # update systematics
	sed -i "/^stats_only/s:^:#:" $atmo
fi

reco_atmo=$root'/../../reconstruction_atmo/*mc.sk4.*.root'
#MC inputs
sed -i "s:^MC_input.*:MC_input\t\"$reco_atmo\":"	$atmo
sed -i "s:^MC_tree_name.*:MC_tree_name\t\"osc_tuple\":"	$atmo

dens=$PWD'/data/PREM_20pts.dat'
prod=$PWD'/data/prod_honda/kam-ally-aa-*.d'
sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$atmo
sed -i "s:production_heights.*:production_heights\t\"$prod\":"	$atmo


#get list of points
#if [ ! -s $global/$pinfo ] ; then
#	echo "Point file does not exist"
#	exit 1
#fi




./bin/oscillation_point $oscc $parm > $root/sensitivity/$pinfo
point=$(cat $root/sensitivity/$pinfo)
point=(${point})

if [ "$SCHED" == "HTCONDOR" ] ; then
	sub=condor_submit
	running="condor_q -autoformat owner jobstatus | grep 2 | wc -l"
	inqueue="condor_q -autoformat owner jobstatus | grep 1 | wc -l"
	echo Launching $NJOBS jobs with HTCondor
elif [ "$SCHED" == "SLURM" ] ; then
	sub=sbatch
	running="squeue -h -r -u $USER -o "%u %t"| grep R  | wc -l"
	inqueue="squeue -h -r -u $USER -o "%u %t"| grep PD | wc -l"
	echo Launching $NJOBS jobs with Slurm
fi

#ready to loop over points
for t in "${point[@]}" ; do

	output=$root/sensitivity/$tname$t
	mkdir -p $output
	sed -i "s:^point.*:point\t$t:" $card
	sed -i "s:^output.*:output\t\"$output/SpaghettiSens.T2HK.root\":" $card

	rm -f $output/L*log
	scriptname=$output/R$nameExec.$t.sub

	this=$output/this_sensitivity.card
	cp $card $this
	sed -i "s:^fit_parameters.*:fit_parameters\t\"$fitc\":" $this
	sed -i "s:^oscillation_parameters.*:oscillation_parameters\t\"$oscc\":" $this
	sed -i "s:^beam_parameters.*:beam_parameters\t\"$beam\":" $this
	sed -i "s:^atmo_parameters.*:atmo_parameters\t\"$atmo\":" $this

	## send as many jobs as files
	if [ "$SCHED" == "HTCONDOR" ] ; then
		cat > $scriptname << EOF
# script submission for condor
# sumbit with --
#	$sub $scriptname

executable		= $Sens
arguments		= \$(Process) $NJOBS $this
getenv			= True
#requirements		= HasAVXInstructions
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $output/L$nameExec.\$(Process).log
error			= $output/L$nameExec.\$(Process).log

queue $NJOBS

EOF
	elif [ "$SCHED" == "SLURM" ] ; then
		cat > $scriptname << EOF
#! /bin/bash
# script submission for SLURM
# sumbit with --
#	$sub $scriptname

#SBATCH --array=0-$((NJOBS - 1))
#SBATCH --job-name=$nameExec
#SBATCH -o $output/L$nameExec.%a.log
#SBATCH -p nms_research,shared
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1

srun $Sens \$SLURM_ARRAY_TASK_ID $NJOBS $this

EOF
	fi

	echo Submitting point $tname$t
	$sub $scriptname

	#wait until jobs are finished before moving to next point
	#but return if last point
	jobs_run=$(eval $running)
	jobs_que=$(eval $inqueue)
	echo $jobs_run and $jobs_que
	if [ $t -ne ${point[${#point[@]} - 1]} ] ; then
		echo not last point.. and running $jobs_run and waiting $jobs_que
		while [ $jobs_run -gt $MAX_JOBS ] || [ $jobs_que -gt $MAX_QUEUE ] ; do
			echo 'waiting 5min...'
			sleep 300
			jobs_run=$(eval $running)
			jobs_que=$(eval $inqueue)
		done
	fi
	sleep 10
done
