#! /bin/bash

#oscillate files from GlobalOsc

Valid=bin/validation
card=/home/tboschi/OscAna/SuperHK/cards/validate.card

root=/home/tboschi/data/
MH=NH

cc=false
while getopts 'r:c' flag; do
	case "${flag}" in
		r) root="${OPTARG}" ;;
		c) cc=true ;;
		*) exit 1 ;;
	esac
done

root=${root%/}

if [ "$cc" = true ]
then
	outdir=$root/validation_corr/
else
	outdir=$root/validation/
fi

mkdir -p $outdir
sed -i "s:base.*:base \"$outdir\":" $card

sE_FHC=$root/systematics/FHC1Re.fij.t2k_spline.root
sE_RHC=$root/systematics/RHC1Re.fij.t2k_spline.root
sM_FHC=$root/systematics/FHC1Rmu.fij.t2k_spline.root
sM_RHC=$root/systematics/RHC1Rmu.fij.t2k_spline.root
matrix=$root/systematics/combinedmatrix.root

sed -i "s:systematic_E_FHC.*:systematic_E_FHC \"$sE_FHC\":" $card
sed -i "s:systematic_E_RHC.*:systematic_E_RHC \"$sE_RHC\":" $card
sed -i "s:systematic_M_FHC.*:systematic_M_FHC \"$sM_FHC\":" $card
sed -i "s:systematic_M_RHC.*:systematic_M_RHC \"$sM_RHC\":" $card

if [ "$cc" = true ]
then
	sed -i "s:correlation.*:correlation \"$matrix\":" $card
else
	sed -i "s:correlation.*:correlation \"\":" $card
fi

$Valid $card #&> /dev/null
