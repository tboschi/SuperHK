#! /bin/bash

PREFIX=$PWD/$1
OA="http://hep.lancs.ac.uk/~tdealtry/oa"

echo "Bulding folder structure in $PREFIX"

mkdir_wget () {
#1 is path, 2 is folder on OA, 3 is file type
	echo Creating $1
	mkdir $PREFIX/$1 -p
	cd $PREFIX/$1

	echo Downloading $3 from $OA/$2
	wget --accept=$3 -r -l1 -nd -np -e robots=off $OA/$2/
	cd -
}

prepare () {
#1 is folder name, 2 and 3 are matrices names
	root=$PREFIX/$1
	shift
	./prepare_systematics.sh -r $root $@
}

mkdir_wget errorstudy/global/reconstruction 190619 'enu_erec_*.root' 

mkdir_wget errorstudy/0/systematics  190829/0 '*.root' 	#ok

#mkdir_wget errorstudy/1a/systematics 190829/1a '*.root' #ok
#mkdir_wget errorstudy/1b/systematics 190829/1b '*.root' #ok
#mkdir_wget errorstudy/2a/systematics 190829/2a '*.root' #ok
#mkdir_wget errorstudy/2b/systematics 190829/2b '*.root' #ok
#
#mkdir_wget errorstudy/6a/systematics 190905/6a '*.root' #ok
#mkdir_wget errorstudy/7a/systematics 190905/7a '*.root' #ok
#mkdir_wget errorstudy/67/systematics 190905/67 '*.root' #ok
#
#mkdir_wget errorstudy/8/systematics  190619/8  '*.root' #ok
#
#mkdir_wget errorstudy/9/systematics  190829/9  '*.root' #ok
#mkdir_wget errorstudy/10/systematics 190829/10 '*.root' #ok
#mkdir_wget errorstudy/NC/systematics 190829/NC '*.root' #ok
#
#mkdir_wget errorstudy/11a/systematics 190619/11a '*.root' #ok
#mkdir_wget errorstudy/11b/systematics 190619/11b '*.root' #ok
#
#mkdir_wget errorstudy/flux_lukas/systematics 190619/flux_lukas '*.root' #ok
#
#mkdir_wget errorstudy/nuenorm_5_corr/systematics 190712/nuenorm '*nue5.t2k.root' #ok
#mkdir_wget errorstudy/nuenorm_5_anti/systematics 190712/nuenorm '*nue5.t2k.root' #ok
#create matrices for nuenorm


banff="DataFit_Postfit_2018_final_v1_180907_sk_eb_valor_order_all_plus_scc.root"
skdetfsi="SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root"

#prepare systematics
prepare 0 $banff $skdetfsi

#prepare 1a "banff_00_nonuenumucc_remove.root" $skdetfsi
#prepare 1b "banff_01_nonuenumucc_remove.root" $skdetfsi
#prepare 2a "banff_02_nonuebarnumubarcc_remove.root" $skdetfsi
#prepare 2a "banff_03_nonuebarnumubarcc_remove.root" $skdetfsi
#
#prepare 6a "banff_nonu2p2h_remove.root" $skdetfsi
#prepare 7a "banff_nonubar2p2h_remove.root" $skdetfsi
#prepare 67 "banff_nonu2p2h_nonubar2p2h_remove.root" $skdetfsi
#
#prepare 8  "banff_increasednueflux.root" $skdetfsi
#
#prepare 9  "banff_06.root" $skdetfsi
#prepare 10 "banff_07.root" $skdetfsi
#prepare NC "banff_14.root" $skdetfsi
#
#prepare 11a $banff "skdetfsi_escale29.root"
#prepare 11b $banff "skdetfsi_escale19.root"
#
#prepare flux_lukas "banff_extraflux" $skdetfsi
#
#prepare nuenorm_5_corr "norm_nue5_nuebar5_corr.root"
#prepare nuenorm_5_anti "norm_nue5_nuebar5_anticorr.root"
