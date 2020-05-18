#! /bin/bash

#oscillate files from GlobalOsc

Nevt=bin/nevents
card=/data/tboschi/HKsens/OscAna/SuperHK/cards/nevents.card

root=/data/tboschi/HKsens/errorstudy
MH=NH

while getopts 'r:g:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		g) global="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

root=${root%/}
outdir=$root/nevents/

global=${global%/}/reconstruction

mkdir -p $outdir
cp $card $outdir/nevents.card
card=$outdir/nevents.card

sed -i "s:base.*:base \"$outdir\":" $card

global_syst=$global'/syst_nuE0_nuE0_FHC.card'
sed -i "s:reco_nuE0_nuE0_FHC.*:reco_nuE0_nuE0_FHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuM0_nuM0_FHC.card'
sed -i "s:reco_nuM0_nuM0_FHC.*:reco_nuM0_nuM0_FHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuM0_nuE0_FHC.card'
sed -i "s:reco_nuM0_nuE0_FHC.*:reco_nuM0_nuE0_FHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuEB_nuEB_FHC.card'
sed -i "s:reco_nuEB_nuEB_FHC.*:reco_nuEB_nuEB_FHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuMB_nuMB_FHC.card'
sed -i "s:reco_nuMB_nuMB_FHC.*:reco_nuMB_nuMB_FHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuMB_nuEB_FHC.card'
sed -i "s:reco_nuMB_nuEB_FHC.*:reco_nuMB_nuEB_FHC\t\"$global_syst\":" $card

global_syst=$global'/syst_nuE0_nuE0_RHC.card'
sed -i "s:reco_nuE0_nuE0_RHC.*:reco_nuE0_nuE0_RHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuM0_nuM0_RHC.card'
sed -i "s:reco_nuM0_nuM0_RHC.*:reco_nuM0_nuM0_RHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuM0_nuE0_RHC.card'
sed -i "s:reco_nuM0_nuE0_RHC.*:reco_nuM0_nuE0_RHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuEB_nuEB_RHC.card'
sed -i "s:reco_nuEB_nuEB_RHC.*:reco_nuEB_nuEB_RHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuMB_nuMB_RHC.card'
sed -i "s:reco_nuMB_nuMB_RHC.*:reco_nuMB_nuMB_RHC\t\"$global_syst\":" $card
global_syst=$global'/syst_nuMB_nuEB_RHC.card'
sed -i "s:reco_nuMB_nuEB_RHC.*:reco_nuMB_nuEB_RHC\t\"$global_syst\":" $card


$Nevt $card #&> /dev/null
