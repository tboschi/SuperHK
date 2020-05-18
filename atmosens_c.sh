#! /bin/bash

#Osc=/home/tboschi/OscAna/hk.atm+beam/submitOsc.sh

Sens=/data/tboschi/HKsens/OscAna/SuperHK/bin/atmofitter
#Sens=/data/tboschi/HKsens/OscAna/SuperHK/special.sh
nameExec=${Sens##*/}
card=/data/tboschi/HKsens/OscAna/SuperHK/cards/atmofit.card
csub=condor_submit

#Stat=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalSpaghettiSensitivityStatOnly
#Fast=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalSpaghettiSensitivityFast
#Fstt=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalSpaghettiSensitivityStatFast

## for fast change point.info
## for stat, comment sys_ lines

root=/data/tboschi
global=/data/tboschi
MH_1=""
MH_2=""
ss=false	#false if systematic fit, true if only stat fit
ff=false	#false if systematic fit, true if fast fit for dCP scan
mtype="correlation"

while getopts '1:2:r:g:m:N:sf' flag; do
	case "${flag}" in
		1) MH_1="${OPTARG}" ;;
		2) MH_2="${OPTARG}" ;;
		r) root="${OPTARG}" ;;
		g) global="${OPTARG}" ;;
		m) mtype="${OPTARG}" ;;
		N) NJOBS="${OPTARG}" ;;
		s) ss=true ;;
		f) ff=true ;;
		*) exit 1 ;;
	esac
done


#global contains specific set, as in global/asim
#root contains specific set, as in root/asim

global=${global%/}

root=${root%/}

#edit card to point to correct systematics
sys_E_FHC=$root/../systematics/atmospheric.root

matrix=$root/../systematics/combinedmatrix.root
mtype="correlation"


#define mass hierarchy to fit
mhfit=$MH_1"_"$MH_2
root=$root/$mhfit

mkdir -p $root/sensitivity/
cp $card $root/sensitivity/sensitivity.card
card=$root/sensitivity/sensitivity.card




##defines input list
input1=$global/$MH_1/global/true.list
input2=$global/$MH_2/global/fit.list
ls $global/$MH_1/global/*.root > $input1
ls $global/$MH_2/global/*.root > $input2

sed -i "s:^true_input.*:true_input\t\"$input1\":" $card
sed -i "s:^fit_input.*:fit_input\t\"$input2\":"   $card





if [ -z $NJOBS ] ; then
	NJOBS=300
fi

echo Launching $NJOBS jobs

#use correct hierrachies
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

#just stats, comment systematics
if [ "$ss" = true ] ; then
	sed -i ":^systematic_:s:^:#:" $card
else
	sed -i ":^#systematic_:s:^#::" $card
fi





if [ "$ff" = true ] ; then
	pinfo="scan.info"
	tname="scan_"
	sed -i "s:^fast.*:fast\t1:" $card
else
	pinfo="point.info"
	tname="point_"
	sed -i "s:^fast.*:fast\t0:" $card
fi

#get list of points
if [ ! -s $global/$pinfo ] ; then
	echo "Point file does not exist"
	exit 1
fi

point=$(cat $global/$pinfo)
point=(${point})

MAX_JOBS=300
MAX_QUEUE=300

#ready to loop over points
for t in "${point[@]}" ; do

	output=$root/sensitivity/$tname$t
	mkdir -p $output
	sed -i "s:^Point.*:Point\t$t:" $card
	sed -i "s:^output.*:output\t\"$output/SpaghettiSens.T2HK.root\":" $card

	rm -f $output/L*log
	scriptname=$output/R$nameExec.$t.sub

	cp $card $output/this_sensitivity.card

	## send as many jobs as files
	cat > $scriptname << EOF
# script submission for condor
# sumbit with --
#	condor_submit $scriptname

executable		= $Sens
arguments		= \$(Process) $NJOBS $output/this_sensitivity.card
getenv			= False
requirements		= HasAVXInstructions
environment		= "LD_LIBRARY_PATH=$LD_LIBRARY_PATH PATH=$PATH HOME=$HOME"
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $output/L$nameExec.\$(Process).log
error			= $output/L$nameExec.\$(Process).log

queue $NJOBS

EOF

	echo Submitting for point $tname$t
	$csub $scriptname

	#wait until jobs are finished before moving to next point
	#but return if last point
	running=$(condor_q -run -format "%s\n" cmd | grep $nameExec | wc -l)
	inqueue=$(condor_q -all -format "%s\n" cmd | grep $nameExec | wc -l)
	inqueue=$(expr $inqueue - $running)
	if [ $t -ne ${point[${#point[@]} - 1]} ] ; then
		echo not last point.. and running $running and waiting $inqueue
		while [ $running -gt $MAX_JOBS -a $inqueue -gt $MAX_QUEUE ] ; do
			echo 'waiting 5min...'
			sleep 300
			running=$(condor_q -run -format "%s\n" cmd | grep $nameExec | wc -l)
			inqueue=$(condor_q -all -format "%s\n" cmd | grep $nameExec | wc -l)
			inqueue=$(expr $inqueue - $running)
		done
	fi
done
