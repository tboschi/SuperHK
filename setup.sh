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
	wget --accept=$3 -nd -np -e robots=off $OA/$2/
	cd -
}

mkdir_wget errorstudy/global/reconstruction 190619 'enu_erec_*.root' 

mkdir_wget errorstudy/0/systematics  190829/0 '*.root' 	#ok

mkdir_wget errorstudy/1a/systematics 190829/1a '*.root' #ok
mkdir_wget errorstudy/1b/systematics 190829/1b '*.root' #ok
mkdir_wget errorstudy/2a/systematics 190829/2a '*.root' #ok
mkdir_wget errorstudy/2b/systematics 190829/2b '*.root' #ok

mkdir_wget errorstudy/6a/systematics 190905/6a '*.root' #ok
mkdir_wget errorstudy/7a/systematics 190905/7a '*.root' #ok
mkdir_wget errorstudy/67/systematics 190905/67 '*.root' #ok

mkdir_wget errorstudy/8/systematics  190619/8  '*.root' #ok

mkdir_wget errorstudy/9/systematics  190829/9  '*.root' #ok
mkdir_wget errorstudy/10/systematics 190829/10 '*.root' #ok
mkdir_wget errorstudy/NC/systematics 190829/NC '*.root' #ok

mkdir_wget errorstudy/11a/systematics 190619/11a '*.root' #ok
mkdir_wget errorstudy/11b/systematics 190619/11b '*.root' #ok

mkdir_wget errorstudy/flux_lukas/systematics 190619/flux_lukas '*.root' #ok

mkdir_wget errorstudy/nuenorm_5_corr/systematics 190712/nuenorm '*nue5.t2k.root' #ok
mkdir_wget errorstudy/nuenorm_5_anti/systematics 190712/nuenorm '*nue5.t2k.root' #ok
#create matrices for nuenorm
