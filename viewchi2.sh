#! /bin/bash

#oscillate files from GlobalOsc

View=/data/tboschi/HKsens/OscAna/SuperHK/bin/viewchi2
card=/data/tboschi/HKsens/OscAna/SuperHK/cards/viewchi2.card

root=/home/tboschi/data/

name=""
while getopts 'r:g:n:' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		g) global="${OPTARG}" ;;
		n) name="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

root=${root%/}
#outdir=$root/runchi2/

input=$global'*.root'
global=${global%/*}/../../../reconstruction

#mkdir -p $outdir
#cp $card $outdir/runchi2.card
#card=$outdir/runchi2.card

sed -i "s:input.*:input \"$input\":" $card
sed -i "s:base.*:base \"$name\":" $card

sE_FHC=$root/systematics/FHC1Re.fij.t2k_spline.root
sE_RHC=$root/systematics/RHC1Re.fij.t2k_spline.root
sM_FHC=$root/systematics/FHC1Rmu.fij.t2k_spline.root
sM_RHC=$root/systematics/RHC1Rmu.fij.t2k_spline.root
matrix=$root/systematics/combinedmatrix.root

sed -i "s:systematic_E_FHC.*:systematic_E_FHC \"$sE_FHC\":" $card
sed -i "s:systematic_E_RHC.*:systematic_E_RHC \"$sE_RHC\":" $card
sed -i "s:systematic_M_FHC.*:systematic_M_FHC \"$sM_FHC\":" $card
sed -i "s:systematic_M_RHC.*:systematic_M_RHC \"$sM_RHC\":" $card
sed -i "s:correlation.*:correlation \"$matrix\":" $card

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


$View $card #&> /dev/null
