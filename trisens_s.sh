#! /bin/bash

usage="$0 -r root_folder -1 [NH | IH] -2 [NH | IH] [-N num_jobs] [-s] 
              [-t fraction] [-f scan_mode] [-m matrix] [-v verbosity] [-h]

Submit fitter jobs to the scheduler.

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
    -h		    Show usage"

Sens=$PWD/cross-fitter.sh
nameExec=${Sens##*/}
nameExec=${nameExec%.*}
sub=sbatch

card=$PWD/cards/multi.card
fitc=$PWD/cards/fit_options.card
oscc=$PWD/cards/oscillation.card
beam=$PWD/cards/beam_sample.card
atmo=$PWD/cards/atmo_sample.card

#Stat=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalSpaghettiSensitivityStatOnly
#Fast=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalSpaghettiSensitivityFast
#Fstt=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalSpaghettiSensitivityStatFast

## for fast change point.info
## for stat, comment sys_ lines

root=$PWD/errorstudy
#global=/data/tboschi
MH_1=""
MH_2=""
ss=false	#false if systematic fit, true if only stat fit
ff=""	#empty if systematic fit, CPV if fast fit for dCP scan
mtype="correlation"
verb="1"

while getopts 'r:1:2:N:t:m:sf:v:h' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		2) MH_2="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
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

#global contains specific set, as in global/asim
#root contains specific set, as in root/asim

#global=${global%/}

# add PWD 
if [[ "$root" != /* ]]
then
	root=$PWD/${root%/}
else
	root=${root%/}
fi


#define mass hierarchy to fit
mhfit=$MH_1"_"$MH_2
root=$root/$mhfit

mkdir -p $root/sensitivity/
cp $card $fitc $oscc $beam $atmo $root/sensitivity/
card=$root/sensitivity/${card##*/}
fitc=$root/sensitivity/${fitc##*/}
oscc=$root/sensitivity/${oscc##*/}
beam=$root/sensitivity/${beam##*/}
atmo=$root/sensitivity/${atmo##*/}

if [ -n $verb ] ; then
	sed -i "s:verbose.*:verbose\t$verb:" $card $fitc $oscc $beam $atmo
fi

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

echo Launching $NJOBS jobs


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

MAX_JOBS=500
MAX_QUEUE=500

#ready to loop over points
for t in "${point[@]}" ; do

	output=$root/sensitivity/$tname$t
	mkdir -p $output
	sed -i "s:^point.*:point\t$t:" $card
	sed -i "s:^output.*:output\t\"$output/SpaghettiSens.T2HK.root\":" $card

	rm -f $output/L*log
	scriptname=$output/R$nameExec.$t.sub

	cp $card $output/this_sensitivity.card
	this=$output/this_sensitivity.card
	sed -i "s:^fit_parameters.*:fit_parameters\t\"$fitc\":" $this
	sed -i "s:^oscillation_parameters.*:oscillation_parameters\t\"$oscc\":" $this
	sed -i "s:^beam_parameters.*:beam_parameters\t\"$beam\":" $this
	sed -i "s:^atmo_parameters.*:atmo_parameters\t\"$atmo\":" $this

	## send as many jobs as files
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
#environment		= "LD_LIBRARY_PATH=$LD_LIBRARY_PATH PATH=$PATH HOME=$HOME"

	echo Submitting for point $tname$t
	$sub $scriptname

	#wait until jobs are finished before moving to next point
	#but return if last point
	running=$(squeue -u $USER -o "%.${#nameExec}j"| grep $nameExec | wc -l)
	#inqueue=$(squeue -u $USER | grep $nameExec | wc -l)
	if [ $t -ne ${point[${#point[@]} - 1]} ] ; then
		echo not last point.. and running $running
		while [ $running -gt $MAX_JOBS ] ; do
			#-a $inqueue -gt $MAX_QUEUE ] ; do
			echo 'waiting 5min...'
			sleep 300
			running=$(squeue -u $USER -o "%.${#nameExec}j"| grep $nameExec | wc -l)
			#inqueue=$(squeue -u $USER | grep $nameExec | wc -l)
		done
	fi
	sleep 10
done
