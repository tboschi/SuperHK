#! /bin/bash

#oscillate files from GlobalOsc

comp=bin/compare
card=/data/tboschi/HKsens/OscAna/SuperHK/cards/compare.card

root=/data/tboschi/HKsens/errorstudy
while getopts 'g:v:' flag; do
	case "${flag}" in
		g) global="${OPTARG}" ;;
		v) valor="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done


outdir=${global%/}/prediction
global=${global%/}/reconstruction
valor=${valor%/}

mkdir -p $outdir
cp $card $outdir/compare.card
card=$outdir/compare.card

sed -i "s:input.*:input \"$valor\":" $card
sed -i "s:base.*:base \"$outdir/\":" $card


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

outfile=$outdir/comparison.tex
rm -f $outfile

echo -e '\documentclass[aspectratio=169]{beamer}' >> $outfile
echo -e '\\usepackage{booktabs}\n' >> $outfile
echo -e '\setbeamertemplate{navigation symbols}{}' >> $outfile
echo -e '\setbeamerfont{frametitle}{size=\\small}' >> $outfile

echo -e >> $outfile
echo -e '\\begin{document}' >> $outfile

dm21="7.53E-5"
dm32="2.509E-3"
t12="0.303977"
for t23 in "0.528" "0.51" "0.6"; do
	for cpv in "-1.601" "0" "-1.570796"; do
		for t13 in "0.0217379" "0.0204168" "0.0230304"; do

			echo using $dm21 $dm32 $t12 $t13 $t23 $cpv
			sed -i "s:M12.*:M12\t\"s$dm21\":" $card
			sed -i "s:M23.*:M23\t\"s$dm32\":" $card
			sed -i "s:S12.*:S12\t\"s$t12\":" $card
			sed -i "s:S13.*:S13\t\"s$t13\":" $card
			sed -i "s:S23.*:S23\t\"s$t23\":" $card
			sed -i "s:dCP.*:dCP\t\"s$cpv\":" $card

			$comp $card > outlog
			output=$(tail -n1 outlog)

			echo -e '\\input{"'$output'"}' >> $outfile
		done
	done
done

echo -e >> $outfile
echo -e '\\end{document}' >> $outfile

cd $outdir

pdflatex $outfile

cd -
