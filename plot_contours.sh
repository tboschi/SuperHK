#!/bin/bash

#############################
#### user defined options ###
#############################


name="one_two_three"

model=("errorstudy/one/NH_NH/contours/point_xxx.root"
       "errorstudy/two/NH_NH/contours/point_yyy.root"
       "errorstudy/three/NH_NH/contours/point_zzz.root")

nickn=("one_xxx"
       "two_yyy"
       "three_zzz")




###################################
#### do not modify script below ###
###################################


Chi2=bin/dropchi2

lims="plot/limits.gpl"
olim="plot/.limits.gpl"
parms=("M12" "M23" "S12" "S13" "S23" "CP")

echo Creating dir plot/$name
mkdir -p plot/$name
card=plot/$name/$name.card

rm $lims $olim

for i in ${!model[*]}; do 
	m="${model[$i]}"
	oscc=${m%/*}
	oscc=${m/contours/sensitivity}/oscillation.card

	echo "pi = 4.*atan(1.)" > $olim
	for p in "${parms[@]}"
	do
		echo "" >> $lims
		cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "x0_" par " = " $2}' >> $olims
		cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "x1_" par " = " $3}' >> $olims
		cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "nn_" par " = " $4}' >> $olims
	done
	if [ ! -f $lims ] ; then # lims file does not exist
		cp $olim $lims
	elif ! diff $lims $olim ; then
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

./bin/plot_bundle plot/$name/$name'_'*.tex
