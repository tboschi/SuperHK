#! /bin/bash

#run GlobalOsc efficiently

Exec=/data/tboschi/HKsens/OscAna/hk.atm+beam/bin/GlobalOsc
csub=condor_submit

global=""

#global already contains hkdr or asim
while getopts 'inf:g:' flag; do
	case "${flag}" in
		i) MH=IH ;;
		n) MH=NH ;;
		f) NFILES="${OPTARG}" ;;
		g) global="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

global=${global%/}

mkdir -p $global/$MH/global

cp $global/*.card $global/$MH/global/global.card
global=$global/$MH/global
card=$global/global.card

#name=$(grep OutputOscFile $card | cut -d'"' -f2)
#output=${name%.*}

#echo sed -i "s:OutputOscFile.*:OutputOscFile \"$output.root\":" $card

if [[ $MH = "NH" ]]; then
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 0:" $card
else
	sed -i "s:InvertedHierarchy.*:InvertedHierarchy 1:" $card
fi

nameExec=${Exec##*Global}
rm -f $global/L$nameExec*log
scriptname=$global/R$nameExec.sub
MAX_JOBS=300

cat > $scriptname << EOF
# script submission for condor
# sumbit with --
#	condor_submit $scriptname


executable		= $Exec
arguments		= \$(Process) $NFILES $global $card
getenv			= True
should_transfer_files	= IF_NEEDED
when_to_transfer_output	= ON_EXIT
initialdir		= $PWD
output			= $global/L$nameExec.\$(Process).log
error			= $global/L$nameExec.\$(Process).log

queue $NFILES
EOF

$csub $scriptname

echo "global finished"
