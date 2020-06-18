#! /bin/bash

usage="$0 [-h] [-f] [-p prefix=\"errorstudy\"]
setup folder structure in working directory
    -h           print this message
    -f           force file overwriting
    -p prefix    specify absolute path for folder structure
    		 default prefix is $PWD/errorstudy"

name=errorstudy
PREFIX=$PWD/errorstudy
force=false
while getopts 'fp:' flag; do
	case "${flag}" in
		f) force=true ;;
		p) PREFIX="${OPTARG}" ;;
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



newbuild=""
PREFIX=${PREFIX%/}

OA="http://hep.lancs.ac.uk/~tdealtry/oa"

echo "Bulding folder structure in $PREFIX"

build () { #1 is path
	if [ "$(ls -A $PREFIX/$1/*.root)" ] && [ "$force" = false ] ; then
		echo Folder $1 not empty.
		echo Call ./setup.sh -f to force overwriting. Now skipping.
		return 1
	fi

	echo Creating $1
	newbuild=$PREFIX/$1
	mkdir $newbuild -p
}

down () { #1 is folder on OA, 2 is file type
	if [ "$(ls -A $newbuild)" ] && [ "$force" = false ] ; then
		return 1
	fi

	mkdir .tmp
	cd .tmp
	echo Downloading $2 from $OA/$1
	echo wget --accept=$2 -r -l1 -nd -np -e robots=off $OA/$1/
	wget --accept=$2 -r -l1 -nd -np -e robots=off $OA/$1/
	mv *.* $newbuild
	cd ..
	rm -rf .tmp
}

create_card () {
#1 is card name, 2 is file name
	card=$PREFIX/global/reconstruction/$1
	file=$PREFIX/global/reconstruction/$2
	cat > $card << EOF
# card file automatically generated

# scale factor, 1 for full statistics, do not touch
scale		1

#TH2D objects to use inside the file
Ring_E_CCQE	"enu_erec_1Re_CCQE"
Ring_E_CCnQE	"enu_erec_1Re_CCnQE"
Ring_E_NC	"enu_erec_1Re_NC"
Ring_M_CCQE	"enu_erec_1Rmu_CCQE"
Ring_M_CCnQE	"enu_erec_1Rmu_CCnQE"
Ring_M_NC	"enu_erec_1Rmu_NC"


reco_path "$file"
EOF
}

prepare () {
#1 is folder name, 2 and 3 are matrices names
	root=$PREFIX/$1
	shift
	./prepare_systematics.sh -r $root $@
}

build global/reconstruction; down 190619 'enu_erec_*.root' 
create_card "syst_nuE0_nuE0_FHC.card" "enu_erec_nue_x_nue.root"
create_card "syst_nuE0_nuE0_RHC.card" "enu_erec_nue_crs_nue_anu.root"
create_card "syst_nuEB_nuEB_FHC.card" "enu_erec_nuebar_x_nuebar.root"
create_card "syst_nuEB_nuEB_RHC.card" "enu_erec_nuebar_crs_nuebar_anu.root"
create_card "syst_nuM0_nuE0_FHC.card" "enu_erec_numu_x_nue.root"
create_card "syst_nuM0_nuE0_RHC.card" "enu_erec_numu_crs_nue_anu.root"
create_card "syst_nuM0_nuM0_FHC.card" "enu_erec_numu_x_numu.root"
create_card "syst_nuM0_nuM0_RHC.card" "enu_erec_numu_crs_numu_anu.root"
create_card "syst_nuMB_nuEB_FHC.card" "enu_erec_numubar_x_nuebar.root"
create_card "syst_nuMB_nuEB_RHC.card" "enu_erec_numubar_crs_nuebar_anu.root"
create_card "syst_nuMB_nuMB_FHC.card" "enu_erec_numubar_x_numubar.root"
create_card "syst_nuMB_nuMB_RHC.card" "enu_erec_numubar_crs_numubar_anu.root"
cp "data/binning.card" $PREFIX/global/reconstruction
cp -r data/asim $PREFIX/global

build 0/systematics; down 190619/0 '*.root' 	#ok

build 1a/systematics; down 190829/1a '*.root' #ok
build 1b/systematics; down 190829/1b '*.root' #ok
build 2a/systematics; down 190829/2a '*.root' #ok
build 2b/systematics; down 190829/2b '*.root' #ok

build 6a/systematics; down 190905/6a '*.root' #ok
build 7a/systematics; down 190905/7a '*.root' #ok
build 67/systematics; down 190905/67 '*.root' #ok

build 8/systematics ; down 190619/8  '*.root' #ok

build 9/systematics ; down 190829/9  '*.root' #ok
build 10/systematics; down 190829/10 '*.root' #ok
build NC/systematics; down 190829/NC '*.root' #ok

build 11a/systematics; down 190619/11a '*.root' #ok
build 11b/systematics; down 190619/11b '*.root' #ok

build flux_lukas/systematics; down 190619/flux_lukas '*.root' #ok

build nuenorm_5_corr/systematics; down 190712/nuenorm '*nue5.t2k.root' #ok
				  down 191028 '*nue5*_corr.root' #ok
build nuenorm_5_anti/systematics; down 190712/nuenorm '*nue5.t2k.root' #ok
				  down 191028 '*nue5*_anticorr.root' #ok


banff="DataFit_Postfit_2018_final_v1_180907_sk_eb_valor_order_all_plus_scc.root"
skdetfsi="SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root"

##prepare systematics
prepare 0 $banff $skdetfsi

prepare 1a "banff_00_nonuenumucc_remove.root" $skdetfsi
prepare 1b "banff_01_nonuenumucc_remove.root" $skdetfsi
prepare 2a "banff_02_nonuebarnumubarcc_remove.root" $skdetfsi
prepare 2b "banff_03_nonuebarnumubarcc_remove.root" $skdetfsi

prepare 6a "banff_nonu2p2h_remove.root" $skdetfsi
prepare 7a "banff_nonubar2p2h_remove.root" $skdetfsi
prepare 67 "banff_nonu2p2h_nonubar2p2h_remove.root" $skdetfsi

prepare 8  "banff_increasednueflux.root" $skdetfsi

prepare 9  "banff_06.root" $skdetfsi
prepare 10 "banff_07.root" $skdetfsi
prepare NC "banff_14.root" $skdetfsi

prepare 11a $banff "skdetfsi_escale29.root"
prepare 11b $banff "skdetfsi_escale19.root"

prepare flux_lukas "banff_extraflux.root" $skdetfsi

prepare nuenorm_5_corr "norm_nue5_nuebar5_corr.root"
prepare nuenorm_5_anti "norm_nue5_nuebar5_anticorr.root"
