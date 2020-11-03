#!/bin/bash



#############################
#### user defined options ###
#############################


# output name unique to combination of sets
# this name will be used to refer to each output file
name="beam_test"


# list of files with contours
# these are the output of the contours.sh script
# and should contain histograms with X2 profiles
model=("errorstudy/one/NH_NH/contours/point_xxx.root"
       "errorstudy/two/NH_NH/contours/point_yyy.root"
       "errorstudy/three/NH_NH/contours/point_zzz.root")


# names for displaying above results
# this will be used by gnuplot to label the curves
# in the plots and correspond to the input files above
nickn=("one_xxx"
       "two_yyy"
       "three_zzz")




###################################
#### do not modify script below ###
###################################


usage="Usage: $0 [-p output] [-h]
			
Render plots out of the contours (see contours.sh) and save them to PDF
In order to use this script, the user should open it with their favourite editor
and modify the first block accordingly. All the outputs will be saved under the 
plot folder, where the plotting scripts are located.
Please refer to the documentation if this is not clear.

  parameters
    -p output    define output PDF file where to save all the files
    -h		 print this message and exit
"

pdf=false
out=""
while getopts 'r:d:1:2:N:t:m:sf:p:xv:h' flag; do
	case "${flag}" in
		p) pdf=true
		   out="${OPTARG}" ;;
		h) echo "$usage" >&2
		   exit 0 ;;
		*) printf "illegal option -%s\n" "$OPTARG" >&2
		   echo "$usage" >&2
		   exit 1 ;;
	esac
done


Chi2=bin/dropchi2

lims="plot/limits.gpl"
olim="plot/.limits.gpl"
parms=("M12" "M23" "S12" "S13" "S23" "CP")

echo Creating dir plot/$name
mkdir -p plot/$name
card=plot/$name/$name.card

rm -f $lims $olim

for i in ${!model[*]}; do 
	m="${model[$i]}"
	oscc=${m%/*}
	oscc=${oscc/contours/sensitivity}/oscillation.card

	echo "pi = 4.*atan(1.)" > $olim
	for p in "${parms[@]}"
	do
		echo "" >> $lims
		cat "$oscc" | awk -v par=$p '/^parm/ && $0~par {gsub('/,/', "", $2); print "x0_" par " = " $2}' >> $olim
		cat "$oscc" | awk -v par=$p '/^parm/ && $0~par {gsub('/,/', "", $3); print "x1_" par " = " $3}' >> $olim
		cat "$oscc" | awk -v par=$p '/^parm/ && $0~par {gsub('/,/', "", $4); print "nn_" par " = " $4}' >> $olim
	done

	cp $olim $lims
	if [ -f $lims ] && ! diff $lims $olim; then # lims file does not exist
		echo oscillation cards are different! Cannot process "$m" with the rest
		exit 1
	fi

	echo "${nickn[$i]}" \"$m\" >> $card
done

echo Extract chi2 information
$Chi2 $card


if [ $name ] ; then
	x2file=($(ls X2*.dat))

	for f in ${x2file[@]}; do
		echo $f
		mv $f plot/$name/$name'_'$f
	done
fi

option="dataset=\"$name\""

plot=( "plot_dCP.gpl"
       "plot_M23.gpl"
       "plot_S23.gpl"
       "plot_S13.gpl"
       "cont_dCP_M23.gpl"
       "cont_S13_dCP.gpl"
       "cont_S23_dCP.gpl"
       "cont_S13_M23.gpl"
       "cont_S23_M23.gpl"
       "cont_S13_S23.gpl" )

cd plot

for p in "${plot[@]}" ; do
	gnuplot -e "$option" $p
done

echo mv -f $name"*.*" $name/
mv -f $name*.* $name/

cd ..

if [ "$pdf" == "true" ] ; then
	if [ -z $out ] ; then
		out="plot_bundle.tex"
	fi

	./bin/plot_bundle $out plot/$name/$name'_'*.tex
fi
