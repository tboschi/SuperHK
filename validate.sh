#! /bin/bash

usage="Usage: $0 -r root_folder -1 [NH | IH] [-p true_point]
              [-2 [NH | IH] ] [-f fit_point] [-t fraction]
	      [-n name] [-x] [-v verbosity] [-h]

Create event samples between at given point; if two points are given, X2 is computed too

  parameters
    -r root_folder  main working folder, must have a \"systematics\" sub folder
    -d data_sample  sample to analyse (beam, data, comb)
    -1 [NH | IH]    mass hierarchy of the observed/true data (NH or IH)

  special parameter
    -x		    use existing cards in output folder if they exist
                    this parameter assumes that root_folder is a full formed path
		    to the sensitivity output location; options -d, -1 and -2 are then ignored

  optional parameters
    -p true	    true point, if not given this is the nominal point
    -2 [NH | IH]    mass hierarchy of the expected/to fit data (NH or IH), default NH
    -f fit	    fit point for comparison, if given X2 is also computed
    -t fraction	    specify fraction of data to fit as value [0,1], default 1 (full data)
    -n name	    give a unique prefix name to outputs, default no prefix
    -v verbosity    specify a verbosity value
    -h		    show usage
"


Vali=$PWD/bin/validation
Oscp=$PWD/bin/oscillation_point

card=$PWD/cards/multi.card
fitc=$PWD/cards/fit_options.card
oscc=$PWD/cards/oscillation.card
beam=$PWD/cards/beam_sample.card
atmo=$PWD/cards/atmo_sample.card


root=""
data=""
#global=/data/tboschi
MH_1=""
MH_2=""
tpoint=""
fpoint=""
exist=false
name=""
verb="1"

while getopts 'r:d:1:2:p:f:t:v:n:xh' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		d) data="${OPTARG}" ;;
		1) MH_1="${OPTARG}" ;;
		2) MH_2="${OPTARG}" ;;
		p) tpoint="${OPTARG}" ;;
		f) fpoint="${OPTARG}" ;;
		t) stats="${OPTARG}" ;;
		v) verb="${OPTARG}" ;;
		n) name="${OPTARG}" ;;
		x) exist=true ;;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done

rm -f .reconstruction_files
rm -f .production_files

# add PWD 
if [[ "$root" != /* ]] ; then
	root=$PWD/${root%/}
else
	root=${root%/}
fi

if [ "$exist" == "true" ] ; then
	card=$root/this_sensitivity.card

	if ! [ -s $card ] ; then
		echo ERROR: file $card does not exist. Resubmit without the -x option
		exit 1
	fi
else # must build folders and card files

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

	if [ -z "$tpoint" ] ; then
		tpoint=$($Oscp $oscc)
		echo True point not defined, so nominal point $tpoint will be used
	else
		echo True point is set to "$tpoint"
	fi

	if [ ! -z "$fpoint" ]; then
		echo Fit point is set to "$fpoint"
	fi


	#define mass hierarchy to fit
	root=$root/validation

	if [ -z "$fpoint" ] ; then
		root=$root/$MH_1'_'$tpoint
	else
		root=$root/$MH_1'_'$tpoint'_'$MH_2'_'$fpoint
	fi

	mkdir -p $root/

	# copy cards to output folder
	cp $card $fitc $oscc $beam $atmo $root/
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

	# uncomment point line
	sed -i "/^#point/s:^#::" $card
	sed -i "s:^point.*:point\t$tpoint:" $card

	if [ -z "$fpoint" ] ; then
		# comment fit point line
		sed -i "/^fit_point/s:^:#:" $card
	else
		# uncomment fit point line
		sed -i "/^#fit_point/s:^#::" $card
		sed -i "s:^fit_point.*:point\t$fpoint:" $card
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
		sed -i "s:stats.*:stats\t$stats:" $beam $atmo
	fi


	# it is best that these remain defined
	#beam systematics

	corr_beam=$root/../systematics/combinedmatrix.root
	corr_atmo=$root/../systematics/atmo_corr.root
	mtype="correlation"

	# need to know correlation matrix to know number of systematics
	sed -i "s:^corr_file.*:corr_file\t\"$corr_beam\":" $beam
	sed -i "s:^corr_name.*:corr_name\t\"$mtype\":"  $beam

	#just stats, comment systematics
	sed -i "/^systematic_/s:^:#:" $beam

	reco_beam=$root'/../../reconstruction_beam/syst_*.card'
	sed -i "s:^reco_input.*:reco_input\t\"$reco_beam\":" $beam

	dens=$PWD'/data/DensityProfileTochibora.dat'
	sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$beam


	#atmo systematics

	sys_atmo=$root/../systematics/atmo_fij.root
	sed -i "s:^systematic_file.*:systematic_file\t\"$sys_atmo\":" $atmo
	sed -i "s:^systematic_tree.*:systematic_tree\t\"sigmatree\":" $atmo
	#Atmo systematics
	sed -i "/^#stats_only/s:^#::" $atmo

	reco_atmo=$root'/../../reconstruction_atmo/*mc.sk4.*.root'
	#MC inputs
	sed -i "s:^MC_input.*:MC_input\t\"$reco_atmo\":"	$atmo
	sed -i "s:^MC_tree_name.*:MC_tree_name\t\"osc_tuple\":"	$atmo

	dens=$PWD'/data/PREM_25pts.dat'
	prod=$PWD'/data/prod_honda/kam-ally-aa-*.d'
	sed -i "s:density_profile.*:density_profile\t\"$dens\":"	$atmo
	sed -i "s:production_heights.*:production_heights\t\"$prod\":"	$atmo


fi

output=$root/"$name"
sed -i "s:^output.*:output\t\"$output\":" $card

$Vali $card
