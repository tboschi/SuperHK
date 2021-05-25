#! /bin/bash

usage="
usage: $0 -r <root> -d <sample> -1 <mh1> -2 <mh2> [-x] [<options>]

Submit fitter jobs to a scheduler, HTCondor or Slurm.
The correct binary run is taken from cross-binary.
The main working folder is defined by <root> which must have a \"systematics\"
subfolder containing systematic input files. The <sample> to fit can be one of
the following: \"beam\", \"data\", or \"comb\" for a combined fit.
The mass hierarchy for the true sample is <mh1> and for the fit sample is <mh2>
and they can both be either \"NH\" for normal hierarchy of \"IH\" for inverted
hierarchy. The output folder will be 

	<root>/<mh1>_<mh2>/sensitivity/
	
The script uses already existing card templates in the \"cards\" folder, modifies
them, and copies them in the output folder.

The special option [-x] will use existing card in the output folder, if they exist.
In this case, the parameters <sample>, <mh1>, and <mh2> are ignored.
This option is useful when resubmitting broken files in combination with the
[-p <list>] option (see below).

Optional parameters
    -N <jobs>       number of jobs to submit to the cluster; the default value
                    is 360. There will be <jobs> output files in the end.
    -s		    does not fit systematic parameters.
    -t <stat>	    specify fraction of data to fit as float value between
    		    0 and 1; the default value is 1, i.e. full data.
    -f <scan>       special scan mode. At the moment only \"CPV\" is available
                    and the fitter will change deltaCP true value along the range
    -p <list>       file with list of true points to fit, useful to use in combin-
                    ation with [-x]. If not specified, nominal points are used.
    -w <dir>        specify a different directory for log files as some file shared
    		    systems redirect output differently
    -m <matrix>	    specify matrix name for beam sample; the default one is
    		    \"correlation\"
    -v <verb>       specify a verbosity value where <verb> is an integer number;
    		    the greater the number, the higher the verbosity.
    -h		    print this message.
"

SCHED=""
if condor_q &> /dev/null ; then
	SCHED="HTCONDOR"
elif squeue &> /dev/null ; then
	SCHED="SLURM"
elif qstat &> /dev/null ; then
    SCHED="SGE"
elif qstat -sql &> /dev/null ; then
    SCHED="SGE"
else
	echo "There is neither HTCondor nor Slurm nor qsub \(SGE or PBS\) on this machine. \
		Iam sorry, I cannot help you"
	exit 1
fi

Bin=$PWD/bin/cross-binary
Sens=fitter
Oscp=bin/oscillation_point

card=$PWD/cards/multi.card
fitc=$PWD/cards/fit_options.card
oscc=$PWD/cards/oscillation.card
beam=$PWD/cards/beam_sample.card
atmo=$PWD/cards/atmo_sample.card

MAX_JOBS=500
MAX_QUEUE=500

root=""
logr=""
data=""
#global=/data/tboschi
MH_1=""
MH_2=""
ss=false	#false if systematic fit, true if only stat fit
ff=""		#empty if systematic fit, CPV if fast fit for dCP scan
list=""
exist=false
NJOBS=360
mtype="correlation"
verb="1"

while getopts 'r:d:1:2:N:t:w:m:sf:p:xv:h' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		2) MH_2="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
		w) logr="${OPTARG}" ;;
		d) data="${OPTARG}" ;;
		m) mtype="${OPTARG}" ;;
		N) NJOBS="${OPTARG}" ;;
		t) stats="${OPTARG}" ;;
		s) ss=true ;;
		f) ff="${OPTARG}" ;;
		p) list="${OPTARG}" ;;
		x) exist=true ;;
		v) verb="${OPTARG}" ;;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done

shift $((OPTIND-1))
if [ $# -gt 0 ]; then
	echo "illegal extra arguments \"$@\"" >&2
	echo "$usage" >&2
	exit 1
fi


# add PWD 
if [[ "$root" != /* ]] ; then
	root=$PWD/${root%/}
else
	root=${root%/}
fi

if [ "$exist" == "true" ] ; then
	card=$root/${card##*/}
	fitc=$root/${fitc##*/}
	oscc=$root/${oscc##*/}
	beam=$root/${beam##*/}
	atmo=$root/${atmo##*/}

	if ! [ -s $card ] ; then
		echo ERROR: file $card does not exist. Resubmit without the -x option
		exit 1
	fi
else # must build folders and card files

	if [ -z $root ] ; then
		echo You must specify a root folder with -r
		exit 1
	fi

	if [ -z $MH_1 ] || [ -z $MH_2 ] ; then
		echo You must define mass hirerachy for true and fitted samples with -1 and -2
		exit 1
	fi

	if [ -z $data ] ; then
		echo You must define a sample to fit with -d
		exit 1
	fi


	#define mass hierarchy to fit
	mhfit=$MH_1"_"$MH_2
	upper=$root
	root=$root/$mhfit/sensitivity

	mkdir -p $root

	if [ -n "$logr" ] ; then
		logr=$logr/$mhfit/sensitivity
		mkdir -p $logr
	fi

	# copy cards to output folder
	cp $card $fitc $oscc $beam $atmo $root
	card=$root/${card##*/}
	fitc=$root/${fitc##*/}
	oscc=$root/${oscc##*/}
	beam=$root/${beam##*/}
	atmo=$root/${atmo##*/}

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

	# decide if normal fit or fast fit for sensitivity
	if [ -n "$ff" ] ; then
		sed -i "/^#scan/s:^#::" $card
		sed -i "s:^scan.*:scan\t$ff:" $card
		sed -i "/^#point/s:^#::" $card
	else
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

	corr_beam=$upper/systematics/matrix.root
	mtype="correlation"

	# need to know correlation matrix to know number of systematics
	sed -i "s:^corr_file.*:corr_file\t\"$corr_beam\":" $beam
	sed -i "s:^corr_name.*:corr_name\t\"$mtype\":"  $beam

	sys_E_FHC=$upper/systematics/FHC1Re.fij.spline.root
	sys_E_RHC=$upper/systematics/RHC1Re.fij.spline.root
	sys_M_FHC=$upper/systematics/FHC1Rmu.fij.spline.root
	sys_M_RHC=$upper/systematics/RHC1Rmu.fij.spline.root
	#just stats, comment systematics
	if [ "$ss" = "true" ] ; then
		echo "statistics only fit"
		sed -i "/^systematic_/s:^:#:" $beam
	else # update systematics
		sed -i "/^#systematic_/s:^#::" $beam
		sed -i "s:^systematic_E_FHC.*:systematic_E_FHC \"$sys_E_FHC\":" $beam
		sed -i "s:^systematic_E_RHC.*:systematic_E_RHC \"$sys_E_RHC\":" $beam
		sed -i "s:^systematic_M_FHC.*:systematic_M_FHC \"$sys_M_FHC\":" $beam
		sed -i "s:^systematic_M_RHC.*:systematic_M_RHC \"$sys_M_RHC\":" $beam
	fi

	reco_beam=$(ls $upper'/../reconstruction_beam/syst_'*.card)
	reco_beam=(${reco_beam})
	reco_beam=(${reco_beam[@]/#/\"})	# creates a nice list of files
	reco_beam=(${reco_beam[@]/%/\",})	# really bad looking , but it works!
	reco_beam="${reco_beam[@]}"
	reco_beam=${reco_beam%,}
	sed -i "s:^reco_input.*:reco_input\t$reco_beam:" $beam

	dens=$PWD'/data/DensityProfileTochibora.dat'
	sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$beam


	#atmo systematics

	sys_atmo=$upper/systematics/atmo_fij.root
	sed -i "s:^systematic_file.*:systematic_file\t\"$sys_atmo\":" $atmo
	sed -i "s:^systematic_tree.*:systematic_tree\t\"sigmatree\":" $atmo
	#Atmo systematics
	if [ "$ss" = "true" ] ; then
		echo "statistics only fit"
		sed -i "/^#stats_only/s:^#::" $atmo
	else # update systematics
		sed -i "/^stats_only/s:^:#:" $atmo
	fi

	reco_atmo=$upper'/../reconstruction_atmo/*.sk4.*.root'
	#MC inputs
	sed -i "s:^MC_input.*:MC_input\t\"$reco_atmo\":"	$atmo
	sed -i "s:^MC_tree_name.*:MC_tree_name\t\"osc_tuple\":"	$atmo

	pre_NH=$upper'/../reconstruction_atmo_fixedSeed_lxplus/pre/NH/atmo.*.root'
	pre_IH=$upper'/../reconstruction_atmo_fixedSeed_lxplus/pre/IH/atmo.*.root'
	#pre computed inputs
	sed -i "s:^pre_input_NH.*:pre_input_NH\t\"$pre_NH\":" $atmo
	sed -i "s:^pre_input_IH.*:pre_input_IH\t\"$pre_IH\":" $atmo
	sed -i "s:^pre_tree_name.*:pre_tree_name\t\"atmoTree\":" $atmo

	dens=$PWD'/data/PREM_25pts.dat'
	prod=$(ls $PWD'/data/prod_honda/'*.d)
	prod=(${prod})
	prod=(${prod[@]/#/\"})	# creates a nice list of files
	prod=(${prod[@]/%/\",})	# really bad looking , but it works!
	prod="${prod[@]}"
	prod=${prod%,}
	sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$atmo
	sed -i "s:honda_production.*:honda_production\t$prod:"	$atmo
fi


# decide if normal fit or fast fit for sensitivity
if grep -q '^scan' $card ; then
	pinfo=$ff"_scan.info"
	tname=$ff"_scan_"
	parm="CP"
else
	pinfo="point.info"
	tname="point_"
	parm=""
fi


if [ -s "$list" ] ; then # user provided list of points
	point=$(cat $list)
else
	$Oscp $oscc $parm > $root/$pinfo
	point=$(cat $root/$pinfo)
fi
point=(${point})

case "$SCHED" in 
	"HTCONDOR")
		sub=condor_submit
		generate=$PWD/src/htcondor.submit
		running="condor_q -autoformat owner jobstatus | grep 2 | wc -l"
		inqueue="condor_q -autoformat owner jobstatus | grep 1 | wc -l"
		;;
	"SLURM")
		sub=sbatch
		generate=$PWD/src/slurm.submit
		running="squeue -h -r -u $USER -o \"%u %t\" | grep R  | wc -l"
		inqueue="squeue -h -r -u $USER -o \"%u %t\" | grep PD | wc -l"
		;;
	"PBS")
		sub=qsub
		generate=$PWD/src/pbs.submit
		running="qstat -t -u $USER | grep \" r \" | wc -l"
		inqueue="qstat -t -u $USER | grep \"qw\" | wc -l"
		;;
	"SGE")
		sub=qsub
		generate=$PWD/src/sge.submit
		running="qstat -t -u $USER | grep \" r \" | wc -l"
		inqueue="qstat -t -u $USER | grep \"qw\" | wc -l"
		;;
esac

#ready to loop over points
for t in "${point[@]}" ; do

	output=$root/$tname$t
	echo $output
	mkdir -p $output
	rm -f $output/*.*
	# changing card in upper folder shows which point is currently being fitted
	sed -i "s:^point.*:point\t$t:" $card
	sed -i "s:^output.*:output\t\"$output/SpaghettiSens.root\":" $card

	scriptname=$output/R$Sens.$t.sh

	this=$output/this_sensitivity.card
	cp $card $this
	sed -i "s:^fit_parameters.*:fit_parameters\t\"$fitc\":" $this
	sed -i "s:^oscillation_parameters.*:oscillation_parameters\t\"$oscc\":" $this
	sed -i "s:^beam_parameters.*:beam_parameters\t\"$beam\":" $this
	sed -i "s:^atmo_parameters.*:atmo_parameters\t\"$atmo\":" $this

	if [ -n "$logr" ] ; then
		outlog=$logr/$tname$t
		mkdir -p $outlog
		echo Log files are stored in $outlog
	else
		outlog=$output
	fi

	# this should work for any script
	$generate $Bin $Sens $NJOBS $this $outlog > $scriptname

	echo Submitting point $tname$t \($NJOBS jobs\) with $SCHED
	$sub $scriptname

	#wait until jobs are finished before moving to next point
	#but return if last point
	jobs_run=$(eval $running)
	jobs_que=$(eval $inqueue)
	if [ $t -ne ${point[${#point[@]} - 1]} ] ; then
		echo not last point.. and running $jobs_run and waiting $jobs_que
		while [ $jobs_run -gt $MAX_JOBS ] || [ $jobs_que -gt $MAX_QUEUE ] ; do
			echo 'waiting 5m...'
			sleep 300
			jobs_run=$(eval $running)
			jobs_que=$(eval $inqueue)
		done
	fi
	sleep 10
done

## send as many jobs as files
#if [ "$SCHED" == "HTCONDOR" ] ; then
#	cat > $scriptname << EOF
## script submission for condor
## sumbit with --
##	$sub $scriptname
#
#executable		= $Sens
#arguments		= fitter \$(Process) $NJOBS $this
#getenv			= True
##requirements		= HasAVXInstructions
#should_transfer_files	= IF_NEEDED
#when_to_transfer_output	= ON_EXIT
#initialdir		= $PWD
#output			= $outlog/L$Sens.\$(Process).log
#error			= $outlog/L$Sens.\$(Process).log
#stream_output		= True
#stream_error		= True
#
#queue $NJOBS
#
#EOF
#elif [ "$SCHED" == "SLURM" ] ; then
#	cat > $scriptname << EOF
##! /bin/bash
## script submission for SLURM
## sumbit with --
##	$sub $scriptname
#
##SBATCH --array=0-$((NJOBS - 1))
##SBATCH --job-name=$Sens
##SBATCH -o $outlog/L$Sens.%a.log
##SBATCH -p nms_research,shared
##SBATCH --time=3-0
##SBATCH --cpus-per-task=1
#
#srun $Sens 
#
#EOF
#else
#	args="\$SGE_TASK_ID $NJOBS $this"
#	outs="$outlog/L$Sens.\$SGE_TASK_ID-1.log"
#	./src/sge.submit $Sens $args $outs $((NJOBS - 1)) > $scriptname
#
#	sub=qsub -t 1-$NJOBS -l sps=1 -o ${outlog}/Logfile$Sens.\$TASK_ID.log -e ${outlog}/Logfile$Sens.\$TASK_ID.log
#	cat > $scriptname << EOF
##! /bin/bash
##script submission for qsub
##submit with qsub -t 0-$NJOBS $scriptname
#
##PBS -N point_$t
##PBS -l walltime=36:00:00
#$Sens fitter \$((SGE_TASK_ID-1)) $NJOBS $this 2>&1 | tee ${outlog}/L$Sens.\$((SGE_TASK_ID-1)).log
#
#EOF
#fi
