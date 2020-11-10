#! /bin/bash

PREFIX=errorstudy
HK="https://pprc.qmul.ac.uk/~tboschi/HK"
usage="$0 [options]
[-f] [-p <prefix>=\"errorstudy\"]
Download systeamtic and reconstruction input files from
database at \"$HK\"
and save them in \"$PREFIX\"

Optional parameters
    -p <prefix>  specify absolute path for folder structure
    		 default prefix is \"$PREFIX\"
    -d <db>      specify a different http database from
    		 where to download input files
    -f           force file overwriting
    -h           print this message
"

PREFIX=$PWD/errorstudy
force=false
while getopts 'hfp:d:' flag; do
	case "${flag}" in
		f) force=true ;;
		d) HK="${OPTARG}" ;;
		p) PREFIX="${OPTARG}" ;;
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

if [[ "$PREFIX" != /* ]]
then
	PREFIX=$PWD/${PREFIX%/}
else
	PREFIX=${PREFIX%/}
fi

echo "Bulding folder structure in $PREFIX"

clone () { #1 is path, #2 is folder on web
	if ls $PREFIX/$1/ >/dev/null 2>&1 && [ "$force" = false ] ; then
		echo Folder $1 not empty.
		echo Call with -f option to force overwriting. Now skipping.
		return 1
	fi

	echo Creating $1
	mkdir $PREFIX/$1 -p

	mkdir .tmp
	cd .tmp
	echo wget from $2
	wget --quiet --accept="*.root" -r -l1 -nd -np -e robots=off $2/
	mv *.* $PREFIX/$1
	cd ..
	rm -rf .tmp
}

card () {
#1 is card name, 2 is file name
	card=$PREFIX/$1/$2
	file=$PREFIX/$1/$3
	cat > $card << EOF
# card file automatically generated

# scale factor, 1 for full statistics, do not touch
scale		1

#TH2D objects to use inside the file
ring_E_CCQE	"enu_erec_1Re_CCQE"
ring_E_CCnQE	"enu_erec_1Re_CCnQE"
ring_E_NC	"enu_erec_1Re_NC"
ring_M_CCQE	"enu_erec_1Rmu_CCQE"
ring_M_CCnQE	"enu_erec_1Rmu_CCnQE"
ring_M_NC	"enu_erec_1Rmu_NC"

reco_path "$file"
EOF
}

combine () {	#1
	out=$PREFIX/$1
	./bin/addmatrix $out/matrix.root $out/*.root
}

# atmo reconstruction consists of only MC files
clone reconstruction_atmo $HK/atmo/reco

# beam reconstruction consists of prediction files and card files
clone reconstruction_beam $HK/beam/reco
card reconstruction_beam "syst_nuE0_nuE0_FHC.card" "enu_erec_nue_x_nue.root"
card reconstruction_beam "syst_nuE0_nuE0_RHC.card" "enu_erec_nue_crs_nue_anu.root"
card reconstruction_beam "syst_nuEB_nuEB_FHC.card" "enu_erec_nuebar_x_nuebar.root"
card reconstruction_beam "syst_nuEB_nuEB_RHC.card" "enu_erec_nuebar_crs_nuebar_anu.root"
card reconstruction_beam "syst_nuM0_nuE0_FHC.card" "enu_erec_numu_x_nue.root"
card reconstruction_beam "syst_nuM0_nuE0_RHC.card" "enu_erec_numu_crs_nue_anu.root"
card reconstruction_beam "syst_nuM0_nuM0_FHC.card" "enu_erec_numu_x_numu.root"
card reconstruction_beam "syst_nuM0_nuM0_RHC.card" "enu_erec_numu_crs_numu_anu.root"
card reconstruction_beam "syst_nuMB_nuEB_FHC.card" "enu_erec_numubar_x_nuebar.root"
card reconstruction_beam "syst_nuMB_nuEB_RHC.card" "enu_erec_numubar_crs_nuebar_anu.root"
card reconstruction_beam "syst_nuMB_nuMB_FHC.card" "enu_erec_numubar_x_numubar.root"
card reconstruction_beam "syst_nuMB_nuMB_RHC.card" "enu_erec_numubar_crs_numubar_anu.root"

# clone main error models
clone systematics_beam/T2K $HK/beam/syst/T2K
combine systematics_beam/T2K	#make correlation matrix
clone systematics_beam/HK  $HK/beam/syst/HK
combine systematics_beam/HK	#make correlation matrix

clone systematics_atmo/SK  $HK/atmo/syst/SK

# clone nuenorm models
clone systematics_beam/NUENORM/anti/1 $HK/beam/syst/NUENORM/anti/1
clone systematics_beam/NUENORM/anti/2 $HK/beam/syst/NUENORM/anti/2
clone systematics_beam/NUENORM/anti/3 $HK/beam/syst/NUENORM/anti/3
clone systematics_beam/NUENORM/anti/4 $HK/beam/syst/NUENORM/anti/4
clone systematics_beam/NUENORM/anti/5 $HK/beam/syst/NUENORM/anti/5

clone systematics_beam/NUENORM/corr/1 $HK/beam/syst/NUENORM/corr/1
clone systematics_beam/NUENORM/corr/2 $HK/beam/syst/NUENORM/corr/2
clone systematics_beam/NUENORM/corr/3 $HK/beam/syst/NUENORM/corr/3
clone systematics_beam/NUENORM/corr/4 $HK/beam/syst/NUENORM/corr/4
clone systematics_beam/NUENORM/corr/5 $HK/beam/syst/NUENORM/corr/5
